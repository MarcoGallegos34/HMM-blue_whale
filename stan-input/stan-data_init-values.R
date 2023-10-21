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
  mu_duration_sigma = c(10,30,30),
  # sigma_duration_alpha = c(3,2.5,15),
  # sigma_duration_beta = c(.022,.008,.03),
  sigma_duration_alpha = c(500,267,845),
  sigma_duration_beta = c(6.25,1.26,6.5),
  
  # surface
  mu_surface_mean = c(70,86,151),
  mu_surface_sigma = c(6,10,10),
  # sigma_surface_alpha = c(1.06,2.45,4.8),
  # sigma_surface_beta = c(.015,.03,.03),
  sigma_surface_alpha = c(361.25,151.25,661.25),
  sigma_surface_beta = c(5.3125,2.75,9.6),
  
  # maxDepth
  mu_maxDepth_mean = c(32,68,170),
  mu_maxDepth_sigma = c(5,10,5),
  # sigma_maxDepth_alpha = c(10,30,40),
  # sigma_maxDepth_beta = c(2,2,2),
  sigma_maxDepth_alpha = c(320,146.7,720),
  sigma_maxDepth_beta = c(13.3,2.25,12),
  
  # step
  mu_step_mean = c(189, 675, 406),
  mu_step_sigma = c(15, 45, 30),
  # sigma_step_alpha = c(134, 305, 287),
  # sigma_step_beta = c(2, 2, 2),
  sigma_step_alpha = c(399, 251.55, 428.56),
  sigma_step_beta = c(2.98, .825, 1.49),
  
  # angle
  # kappa_alpha = c(1,3.1,.8),
  # kappa_beta = c(1,1,1),
  kappa_alpha = c(125,133.47,80),
  kappa_beta = c(125,43.05,100),
  
  # heading variance
  # a_alpha = c(1,.5,1.7),
  # a_beta = c(1,1,1),
  a_alpha = c(125,125,361.25),
  a_beta = c(125,250,212.5),
  b_alpha = c(245,5.4,50.45),
  b_beta = c(116.6,1,9.34),
  
  # lunges
  # lambda_alpha = c(.7,.005,3.4),
  # lambda_beta = c(1,1,1)
  lambda_alpha = c(245,.0138,1445),
  lambda_beta = c(350,2.77,425)
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
