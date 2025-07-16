# Load required packages
library("INLA")
library("raster")
library("INLAspacetime")
library("sf")

# Set seed for reproducibility
set.seed(201805)

# Load spatial files used to build the mesh and define the domain
greenhouses <- shapefile("/path/to/your/domain_dense_greenhouses.shp")
barriers <- shapefile("/path/to/your/dense_greenhouses.shp")
barriers_sf <- st_as_sf(barriers)  # Convert to sf object for geometry operations
points <- shapefile("/path/to/your/points_clustering_dense_greenhouses.shp")  # Training points
testingpoints <- shapefile("/path/to/your/testingpoints_clustering_dense_greenhouses.shp")  # Testing points

# Loop over simulation settings (currently only one value for 'range')
for (j in c(100)) {
  # Set mesh and barrier parameters
  range <- j  # Spatial range in meters
  bound.outer <- 200  # Outer boundary extension for the mesh
  max.edge <- (1/10) * range  # Maximum edge length for mesh triangles
  
  # Build spatial mesh around training points and greenhouse domain
  mesh = inla.mesh.2d(
    boundary = greenhouses,
    loc = cbind(points$x, points$y),
    max.edge = c(2, 4) * max.edge,
    cutoff = 2,  # Minimum allowed distance between points (2 m for scarce, 9 m for dense, 7 m for moderate)
    crs = CRS("+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs"),
    offset = c(max.edge, bound.outer)
  )
  
  # Identify mesh triangles outside greenhouses (i.e., barriers)
  water.tri <- fm_contains(greenhouses, y = mesh, type = "centroid", ignore.CRS = TRUE)
  num.tri <- length(mesh$graph$tv[, 1])
  barrier.tri <- setdiff(1:num.tri, water.tri[[1]])
  poly.barrier <- inla.barrier.polygon(mesh, barrier.triangles = barrier.tri)
  crs(poly.barrier) <- CRS("+proj=utm +zone=30 +datum=WGS84 +units=m +no_defs")
  
  # Convert mesh nodes to spatial points
  mesh_coords <- mesh$loc[, 1:2]
  mesh_pts <- st_as_sf(data.frame(x = mesh_coords[,1], y = mesh_coords[,2]), coords = c("x", "y"), crs = 32630)
  
  # Extract training point coordinates
  coords <- as.data.frame(matrix(c(points$x, points$y), ncol=2))
  colnames(coords) <- c("x", "y")
  
  # Calculate distance from each mesh node to the nearest barrier edge
  barriers_edges <- st_boundary(barriers_sf)
  mesh_dist_matrix_edges <- st_distance(mesh_pts, barriers_edges)
  mesh_min_dist_to_edge <- apply(mesh_dist_matrix_edges, 1, min)
  
  # Define precision matrix for the barrier SPDE model
  prec <- 1
  barrier.model <- inla.barrier.pcmatern(
    mesh, 
    prior.range = c(range, 0.5), 
    prior.sigma = c(prec, 0.5),
    barrier.triangles = barrier.tri,
    range.fraction = 0.8, 
    enable.INLAspacetime = FALSE
  )
  Q <- inla.rgeneric.q(barrier.model, "Q", theta = c(log(prec), log(range)))
  
  # Sample a spatial field realization from the precision matrix
  set.inla.seed <- 2000
  u <- inla.qsample(n = 1, Q = Q, seed = set.inla.seed)
  u <- u[, 1] + exp(-0.001 * (mesh_min_dist_to_edge)^2)  # Add decay near barriers
  
  ######## Simulate data and evaluate model performance ########
  
  # Data frame to store performance metrics
  sensib <- data.frame(size = numeric(70), MAE = numeric(70), ppvalue = numeric(70))
  iter <- rep(c(10), 1, each = 10)
  ii <- 0
  
  for (i in iter) {
    ii <- ii + 1
    
    # Define full training and testing coordinate matrices
    loc.data <- matrix(c(points$x, points$y), 40, 2)
    test_loc.data <- matrix(c(testingpoints$x, testingpoints$y), 40, 2)
    
    # Subsample training and testing data
    sub_sampleSize <- i
    data <- loc.data[sample(1:nrow(loc.data), sub_sampleSize, replace = FALSE), ]
    test_loc.data <- test_loc.data[sample(1:nrow(test_loc.data), round((i * 25) / 100), replace = FALSE), ]
    
    # Create projection matrices (A) from mesh to data locations
    A.data <- inla.spde.make.A(mesh, as.matrix(data))
    A.test_data <- inla.spde.make.A(mesh, test_loc.data)
    
    # Project spatial field values onto training locations
    u.data <- A.data %*% u
    
    # Simulate observed values with Gaussian noise
    df <- data.frame(data)
    names(df) <- c('locx', 'locy')
    sigma.u <- 3         # Size of spatial effect
    sigma.epsilon <- 1   # Size of Gaussian noise
    df$y <- drop(sigma.u * u.data + sigma.epsilon * rnorm(nrow(df)))
    
    # Build INLA stacks for training and prediction
    stk <- inla.stack(
      data = list(y = df$y), A = list(A.data, 1),
      effects = list(s = 1:mesh$n, intercept = rep(1, nrow(df))),
      remove.unused = FALSE, tag = 'est'
    )
    test_stk <- inla.stack(
      data = list(y = rep(NA, nrow(test_loc.data))), A = list(A.test_data, 1),
      effects = list(s = 1:mesh$n, intercept = rep(1, nrow(test_loc.data))),
      remove.unused = FALSE, tag = 'pred'
    )
    join.stack <- inla.stack(stk, test_stk)  # Combine stacks
    
    # Define barrier SPDE model formula
    formula2 <- y ~ 0 + intercept + f(s, model = barrier.model)
    
    # Fit the model using INLA
    res.barrier <- inla(
      formula2, data = inla.stack.data(join.stack),
      control.predictor = list(compute = TRUE, A = inla.stack.A(join.stack)),
      family = 'gaussian', num.threads = 1,
      control.family = list(hyper = list(prec = list(
        prior = "gaussian", fixed = FALSE, control.inla = list(max.iter = 50),
        param = c(0, 1e6)))),
      control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE, return.marginals.predictor = TRUE),
      control.mode = list(restart = TRUE, theta = c(log(1), log(3), log(range))),
      verbose = FALSE
    )
    
    # Generate test values and assess prediction quality
    u.testdata <- A.test_data %*% u
    test_df <- data.frame(test_loc.data)
    names(test_df) <- c('locx', 'locy')
    sigma.u <- 3
    sigma.epsilon <- 1
    test_df$y <- drop(sigma.u * u.testdata + sigma.epsilon * rnorm(nrow(test_df)))
    
    # Compute posterior predictive p-values
    predicted.p.value <- c()
    for (k in 1:nrow(test_df)) {
      predicted.p.value[k] <- inla.pmarginal(
        q = test_df$y[k],
        marginal = res.barrier$marginals.fitted.values[[sub_sampleSize + k]]
      )
    }
    
    # Store results: MAE and median p-value
    sensib$size[ii] <- i
    sensib$MAE[ii] <- sum(abs(res.barrier$summary.fitted.values$mean[(i + 1):(i + round((i * 25) / 100))] - test_df$y)) / round((i * 25) / 100)
    sensib$ppvalue[ii] <- median(predicted.p.value)
    
    print(ii)
  }
  
  # Save model performance metrics
  write.csv(sensib, file = paste0("/path/to/your/results/output", "range_", j, ".csv"), row.names = FALSE)
}
