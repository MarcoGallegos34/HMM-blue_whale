stan_data = list(
  K=1,
  N =3,
  n = nrow(whales_seq),
  n_ind = nrow(whales_init),
  n_ind_total = 37,
  ID_init = whales_init$ID,
  ID = whales_seq$ID,
  x_duration_init=as.numeric(whales_init$DIVE.TIME),
  x_surface_init=as.numeric(whales_init$SURFACE.TIME),
  x_maxDepth_init = as.numeric(whales_init$MAX.DEPTH),
  x_lunges_init = as.numeric(whales_init$LUNGES),
  x_step_init = as.numeric(whales_init$steps),
  x_angle_init = as.numeric(whales_init$turns),
  x_headVar_init = as.numeric(whales_init$var.head),
  x_exposure_init = as.numeric(whales_init$EXPOSURE),
  x_duration = as.numeric(whales_seq$DIVE.TIME),
  x_surface = as.numeric(whales_seq$SURFACE.TIME),
  x_maxDepth = as.numeric(whales_seq$MAX.DEPTH),
  x_lunges = whales_seq$LUNGES,
  x_step = whales_seq$steps,
  x_angle = whales_seq$turns,
  x_headVar = whales_seq$var.head,
  x_exposure = whales_seq$EXPOSURE,
  # hyperparameters to provid
  
  # duration
  mu_duration_mean = c(140,334,516),
  mu_duration_sigma = c(100,300,300),
  # sigma_duration_alpha = c(3,2.5,15),
  # sigma_duration_beta = c(.022,.008,.03),
  sigma_duration_alpha = c(20,20,33.8),
  sigma_duration_beta = c(.25,0.09433962,.26),
  
  # surface
  mu_surface_mean = c(70,86,151),
  mu_surface_sigma = c(60,50,100),
  # sigma_surface_alpha = c(1.06,2.45,4.8),
  # sigma_surface_beta = c(.015,.03,.03),
  sigma_surface_alpha = c(9.248,67.222,9.9146189),
  sigma_surface_beta = c(.136,1.222,0.1436901),
  
  # maxDepth
  mu_maxDepth_mean = c(32,68,170),
  mu_maxDepth_sigma = c(25,50,100),
  # sigma_maxDepth_alpha = c(10,30,40),
  # sigma_maxDepth_beta = c(2,2,2),
  sigma_maxDepth_alpha = c(180,845,45),
  sigma_maxDepth_beta = c(7.5,13,.75),
  
  # step
  mu_step_mean = c(189, 675, 406),
  mu_step_sigma = c(150, 450, 300),
  # sigma_step_alpha = c(134, 305, 287),
  # sigma_step_beta = c(2, 2, 2),
  sigma_step_alpha = c(8.978, 46.5125, 41.1845),
  sigma_step_beta = c(.067, .1525, .1435),
  
  # angle
  # kappa_alpha = c(1,3.1,.8),
  # kappa_beta = c(1,1,1),
  kappa_alpha = c(.3125,4805,.20),
  kappa_beta = c(.3125,1550,.25),
  
  # heading variance
  # a_alpha = c(1,.5,1.7),
  # a_beta = c(1,1,1),
  a_alpha = c(1.25,.2,14.45),
  a_beta = c(1.25,.4,8.5),
  # b_alpha = c(5.5125,25.3125,12.8),
  # b_beta = c(2.6250,4.6875,8),
  b_alpha = c(5.5125,145.8,12.8),
  b_beta = c(2.6250,27,8),
  
  # lunges
  # lambda_alpha = c(.7,.005,3.4),
  # lambda_beta = c(1,1,1)
  lambda_alpha = c(.6125,.000125,14.45),
  lambda_beta = c(.8750,.025,4.25)
)



init_list_m2 = list(
  # duration
  mu_duration = c(140.0,334.0,516.0),
  #mu_duration = c(100.0,100.0,100.0),
  #log_mu_duration = log(c(140.0,334.0,516.0)),
  log_sigma_duration = log(c(80.0,212.0,130.0)),
  
  # surface
  mu_surface = c(70.0,86.0,151.0),
  # log_mu_surface = log(c(70.0,86.0,151.0)),
  log_sigma_surface = log(c(68.0,55.0,69.0)),
  
  #max depth
  mu_maxDepth = c(32.0,68.0,170.0),
  # log_mu_maxDepth = log(c(32.0,68.0,170.0)),
  log_sigma_maxDepth = log(c(24.0,65.0,60.0)),
  
  #step length
  mu_step = c(189.0,675.0,406.0),
  # log_mu_step = log(c(189.0,675.0,406.0)),
  log_sigma_step = log(c(134.0,305.0,287.0)),
  
  #turning angle
  log_kappa = log(c(1.0,3.1,.8)),
  
  #heading variance
  log_a = log(c(1.0,.5,1.7)),
  log_b = log(c(2.1,5.4,1.6)),
  
  #number of lunges
  log_lambda = log(c(.7,.05,3.4)),
  
  #initial distribution
  init_raw = qlogis(c(1/3,1/3)),
  # init = c(1/3,1/3,1/3),
  # init_raw = matrix(qlogis(1/3),3,2),
  
  # pi_raw = qlogis(c(1/3,1/3)),
  # pi_raw = qlogis(c(.9999)),
  
  #rows tpm
  # tpm1 = matrix(0,3,2),
  # tpm2 = matrix(0,3,2),
  # tpm3 = matrix(0,3,2))
  
  tpm1 = c(1/3,1/3,1/3),
  tpm2 = c(1/3,1/3,1/3),
  tpm3 = c(1/3,1/3,1/3))
# tpm1 = matrix(c(0.0,0.0),1,2),
# tpm2 = matrix(c(0.0,0.0),1,2),
# tpm3 = matrix(c(0.0,0.0),1,2))

init_list = list(
  # duration
  mu_duration = c(140.0,334.0,516.0),
  #mu_duration = c(100.0,100.0,100.0),
  #log_mu_duration = log(c(140.0,334.0,516.0)),
  log_sigma_duration = log(c(80.0,212.0,130.0)),
  
  # surface
  mu_surface = c(70.0,86.0,151.0),
  # log_mu_surface = log(c(70.0,86.0,151.0)),
  log_sigma_surface = log(c(68.0,55.0,69.0)),
  
  #max depth
  mu_maxDepth = c(32.0,68.0,170.0),
  # log_mu_maxDepth = log(c(32.0,68.0,170.0)),
  log_sigma_maxDepth = log(c(24.0,65.0,60.0)),
  
  #step length
  mu_step = c(189.0,675.0,406.0),
  # log_mu_step = log(c(189.0,675.0,406.0)),
  log_sigma_step = log(c(134.0,305.0,287.0)),
  
  #turning angle
  log_kappa = log(c(1.0,3.1,.8)), 
  
  #heading variance
  log_a = log(c(1.0,.5,1.7)),
  log_b = log(c(2.1,5.4,1.6)),
  
  #number of lunges
  log_lambda = log(c(.7,.05,3.4)),
  
  #initial distribution
  # init_raw = qlogis(c(1/3,1/3)),
  init = c(1/3,1/3,1/3),
  
  pi_ = c(1/3,1/3,1/3),
  
  #rows tpm
  # tpm1 = c(1/3,1/3,1/3),
  # tpm2 = c(1/3,1/3,1/3),
  # tpm3 = c(1/3,1/3,1/3))
  tpm1 = c(0.0,0.0),
  tpm2 = c(0.0,0.0),
  tpm3 = c(0.0,0.0))
