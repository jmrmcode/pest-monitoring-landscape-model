# Load required packages
library("INLA")
library("raster")
library("INLAspacetime")
library("sf")

set.seed(201805)

# Load required spatial files for building the mesh
# Replace the paths below with the location of your shapefiles
greenhouses <- shapefile("/path/to/your/domain_dense_greenhouses.shp")  # Greenhouse boundary shapefile
barriers <- shapefile("/path/to/your/dense_greenhouses.shp")            # Internal barriers shapefile
barriers_sf <- st_as_sf(barriers)  # Convert to sf object for easier geometry handling

# Extract a random sample of 40 points outside greenhouses
points <- spsample(greenhouses, n = 40, type = "random")


# Loop over simulation settings (currently only one value for 'range')
nrange <- c(100)
for (j in nrange) {
  # Set the length of the boundary extension
  range <- j # in meters
  bound.outer = 200
  max.edge <- (1/10) * range
  
  # Build the mesh; lengths in meters
  mesh = inla.mesh.2d(boundary = greenhouses,
                      loc = cbind(points$x, points$y),
                      max.edge = c(2, 4) * max.edge,
                      cutoff = 2,
                      crs = CRS("+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs"),
                      offset = c(max.edge, bound.outer))
  
  # Select the triangles of the mesh that are inside greenhouses
  water.tri = fm_contains(greenhouses, y = mesh, type = "centroid", ignore.CRS = TRUE)
  num.tri = length(mesh$graph$tv[, 1])
  barrier.tri = setdiff(1:num.tri, water.tri[[1]])
  poly.barrier = inla.barrier.polygon(mesh, barrier.triangles = barrier.tri)
  crs(poly.barrier) <- CRS("+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs")
  
  # Get mesh coordinates and convert to sf
  mesh_coords <- mesh$loc[, 1:2]
  mesh_pts <- st_as_sf(data.frame(x = mesh_coords[,1], y = mesh_coords[,2]), coords = c("x", "y"), crs = 32630)
  
  # Get training point coordinates
  coords <- as.data.frame(matrix(c(points$x, points$y), ncol=2))
  colnames(coords) <- c("x", "y")
  
  # Extract the polygon boundaries (linestrings)
  barriers_edges <- st_boundary(barriers_sf)
  
  # Distance matrix: each point to each polygon edge
  mesh_dist_matrix_edges <- st_distance(mesh_pts, barriers_edges)
  mesh_min_dist_to_edge <- apply(mesh_dist_matrix_edges, 1, min)
  
  # Define the precision matrix for the Barrier model considering unit marginal variance
  prec <- 1
  barrier.model <- inla.barrier.pcmatern(mesh, prior.range = c(range, 0.5), prior.sigma = c(prec, 0.5),
                                         barrier.triangles = barrier.tri, range.fraction = 0.8, enable.INLAspacetime = FALSE)
  Q <- inla.rgeneric.q(barrier.model, "Q", theta = c(log(prec), log(range)))
  
  # Sample from the precision matrix of the spatial field
  set.inla.seed = 2000
  u = inla.qsample(n = 1, Q = Q, seed = set.inla.seed)
  u = u[, 1] + exp(-0.001 * (mesh_min_dist_to_edge)^2)
  
  ####### Extract a random sample from u using points
  niter <- 10
  sensib <- data.frame(size = numeric(7*niter), MAE = numeric(7*niter), ppvalue = numeric(7*niter))
  iter <- rep(c(niter), 1, each = niter)
  ii = 0
  
  for (i in iter) {
    ii <- ii + 1
    test_points <- spsample(greenhouses, n = round((i*25)/100), type = "random")  # Independent random test sample
    
    loc.data <- matrix(c(points$x, points$y), 40, 2)
    test_loc.data <- matrix(c(test_points$x, test_points$y), round((i*25)/100), 2)
    
    sub_sampleSize <- i
    data <- loc.data[sample(1:nrow(loc.data), sub_sampleSize, replace = FALSE), ]
    
    # test_loc.data <- matrix(c(testingpoints$x, testingpoints$y), 80, 2)  # For cluster sampling
    # test_loc.data <- test_loc.data[sample(1:nrow(test_loc.data), round((i*25)/100), replace = FALSE), ]
    
    A.data <- inla.spde.make.A(mesh, as.matrix(data))
    A.test_data <- inla.spde.make.A(mesh, test_loc.data)
    
    u.data <- A.data %*% u
    
    df <- data.frame(data)
    names(df) <- c("locx", "locy")
    sigma.u <- 3
    sigma.epsilon <- 1
    df$y <- drop(sigma.u * u.data + sigma.epsilon * rnorm(nrow(df)))
    
    stk <- inla.stack(data = list(y = df$y), A = list(A.data, 1),
                      effects = list(s = 1:mesh$n, intercept = rep(1, nrow(df))),
                      remove.unused = FALSE, tag = "est")
    
    test_stk <- inla.stack(data = list(y = rep(NA, nrow(test_loc.data))), A = list(A.test_data, 1),
                           effects = list(s = 1:mesh$n, intercept = rep(1, nrow(test_loc.data))),
                           remove.unused = FALSE, tag = "pred")
    
    join.stack <- inla.stack(stk, test_stk)
    
    formula2 <- y ~ 0 + intercept + f(s, model = barrier.model)
    
    res.barrier <- inla(formula2, data = inla.stack.data(join.stack),
                        control.predictor = list(compute = TRUE, A = inla.stack.A(join.stack)),
                        family = "gaussian", num.threads = 1,
                        control.family = list(hyper = list(prec = list(
                          prior = "gaussian", fixed = FALSE, control.inla = list(max.iter = 50),
                          param = c(0, 1e6)))),
                        control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, return.marginals.predictor = TRUE),
                        control.mode = list(restart = TRUE, theta = c(log(1), log(3), log(range))),
                        verbose = FALSE)
    
    ### Model performance using test data
    u.testdata <- A.test_data %*% u
    test_df <- data.frame(test_loc.data)
    names(test_df) <- c("locx", "locy")
    test_df$y <- drop(sigma.u * u.testdata + sigma.epsilon * rnorm(nrow(test_df)))
    
    predicted.p.value <- numeric(nrow(test_df))
    for (k in 1:nrow(test_df)) {
      predicted.p.value[k] <- inla.pmarginal(q = test_df$y[k], marginal = res.barrier$marginals.fitted.values[[sub_sampleSize + k]])
    }
    
    sensib$size[ii] <- i
    sensib$MAE[ii] <- mean(abs(res.barrier$summary.fitted.values$mean[(i+1):(i+round((i*25)/100))] - test_df$y))
    sensib$ppvalue[ii] <- median(predicted.p.value)
    
    print(ii)
  }
  
  # Save results — update this path as needed
  write.csv(sensib, file = paste0("/path/to/your/results/output", "range_", j, ".csv"), row.names = FALSE)
}
