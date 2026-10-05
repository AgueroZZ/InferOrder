# Exact generating functions, reparameterized by their two-feature arc length.
source('experiments/estimate_intrinsic_m_smooth_v032/auto_m_v034/common.R')
suppressPackageStartupMessages(library(ggplot2))
base <- file.path(study_dir,'exploratory_m5_ordering_b')
out <- file.path(base,'arc_length_visualization')
d <- load_fixed_dataset(manifest[manifest$id=='main_M5_S4_r001',])
# Deterministic illustrative choice, independent of fitted method performance.
features <- head(sort(d$coefficients$feature[d$coefficients$ordering=='B' & !d$coefficients$is_anchor]),2)
stopifnot(identical(features,c('V13','V14')))
coef <- d$coefficients[match(features,d$coefficients$feature),]
evaluate <- function(t,derivative=0L) {
 k <- 2:4; angle <- outer(t,pi*k)
 vapply(seq_len(nrow(coef)),function(j) {
  a <- as.numeric(coef[j,paste0('sine_k',k)])/k^2
  b <- as.numeric(coef[j,paste0('cosine_k',k)])/k^2
  value <- if(derivative==0L) as.numeric(sin(angle)%*%a+cos(angle)%*%b) else
   as.numeric(cos(angle)%*%(a*pi*k)-sin(angle)%*%(b*pi*k))
  if(derivative==0L)value <- value-coef$raw_center[j]
  value/coef$raw_scale[j]
 },numeric(length(t)))
}
t <- seq(0,1,length.out=20001L)
f <- evaluate(t); derivative <- evaluate(t,1L);speed <- sqrt(rowSums(derivative^2))
stopifnot(min(speed)>0)
v <- c(0,cumsum(diff(t)*(head(speed,-1)+tail(speed,-1))/2));L <- tail(v,1)
v_grid <- seq(0,L,length.out=20001L)
t_inverse <- approx(v,t,xout=v_grid)$y
g <- evaluate(t_inverse)
g_derivative <- evaluate(t_inverse,1L)/sqrt(rowSums(evaluate(t_inverse,1L)^2))
truth_t <- d$latent_positions[,match('B',d$ordering_labels)]
reconstruction_error <- max(abs(evaluate(truth_t)-d$signal[,match(features,colnames(d$signal))]))
stopifnot(reconstruction_error<1e-12,max(abs(sqrt(rowSums(g_derivative^2))-1))<1e-12)
# Check the quadrature against a coarser grid and the inverse map against f(t).
coarse <- seq(1,length(t),by=2L)
coarse_L <- sum(diff(t[coarse])*(head(speed[coarse],-1)+tail(speed[coarse],-1))/2)
roundtrip_error <- max(abs(evaluate(approx(v_grid,t_inverse,xout=v)$y)-f))
stopifnot(abs(coarse_L-L)<1e-5,roundtrip_error<1e-5)
write.csv(data.frame(t=t,v=v,v_normalized=v/L,f_V13=f[,1],f_V14=f[,2],speed=speed),file.path(out,'true_t_and_arc_length.csv'),row.names=FALSE)
write.csv(data.frame(v=v_grid,v_normalized=v_grid/L,t=t_inverse,g_V13=g[,1],g_V14=g[,2]),file.path(out,'reparameterized_functions.csv'),row.names=FALSE)
checks <- data.frame(feature_1=features[1],feature_2=features[2],total_length=L,
 min_original_speed=min(speed),max_original_speed=max(speed),
 true_signal_reconstruction_error=reconstruction_error,length_grid_error=abs(coarse_L-L),roundtrip_error=roundtrip_error,
 unit_speed_error=max(abs(sqrt(rowSums(g_derivative^2))-1)))
write.csv(checks,file.path(out,'verification.csv'),row.names=FALSE)
saveRDS(list(coefficients=coef,input_hash=d$input_hash,features=features,selection='First two lexicographically sorted non-anchor B feature names',
 dataset_id=d$row$id,domain=c(0,1),length=L,grid_size=length(t),quadrature='Trapezoidal rule on analytic speed',
 metric='Euclidean distance in the two selected, originally standardized feature coordinates',
 script_sha256=digest::digest(file=file.path(out,'render.R'),algo='sha256')),file.path(out,'provenance.rds'))
colors <- c(V13='#0072B2',V14='#D55E00')
theme_set(theme_minimal(base_size=12))
line_theme <- theme(legend.position='bottom',panel.grid.minor=element_blank(),plot.title=element_text(face='bold',size=13))
make_function_plot <- function(x,y,label,title,subtitle) {
 frame <- data.frame(x=rep(x,2),value=c(y[,1],y[,2]),feature=rep(features,each=length(x)))
 ggplot(frame,aes(x,value,color=feature,linetype=feature))+geom_line(linewidth=.85)+
 scale_color_manual(values=colors)+scale_linetype_manual(values=c(V13='solid',V14='dashed'))+
 coord_cartesian(ylim=range(f))+labs(x=label,y='True feature value',title=title,subtitle=subtitle,color=NULL,linetype=NULL)+line_theme
}
p1 <- make_function_plot(t,f,'Original latent position t','A. True functions on the t scale','Same generating functions as the saved B benchmark')
p2 <- make_function_plot(v_grid,g,'Arc length v','B. Reparameterized functions on the v scale',sprintf('g(v) = f(t(v)); total length L = %.3f; ||g\'(v)|| = 1',L))
mark_t <- seq(0,1,length.out=11)
mark_v <- seq(0,L,length.out=11)
points_t <- evaluate(mark_t)
points_v <- evaluate(approx(v,t,xout=mark_v)$y)
curve <- data.frame(x=f[,1],y=f[,2])
geometry <- function(points,title,subtitle,color,shape) {
 dots <- data.frame(x=points[,1],y=points[,2])
 ggplot(curve,aes(x,y))+geom_path(color='#555555',linewidth=.65)+
 geom_point(data=dots,color=color,shape=shape,size=2.8)+
 geom_point(data=dots[c(1,11),],shape=21,fill='white',color=color,size=4,stroke=1.1)+
 coord_equal(xlim=range(f[,1]),ylim=range(f[,2]))+
 labs(x='V13',y='V14',title=title,subtitle=subtitle)+line_theme
}
p3 <- geometry(points_t,'C. The curve with equally spaced t markers','11 points: t = 0, 0.1, ..., 1','#0072B2',16)
p4 <- geometry(points_v,'D. The SAME curve with equal arc-length markers','11 points: v/L = 0, 0.1, ..., 1','#D55E00',17)
p5 <- ggplot(data.frame(t=t,v=v),aes(t,v))+geom_abline(intercept=0,slope=L,color='#999999',linetype=2)+
 geom_line(color='#009E73',linewidth=.9)+
 labs(x='Original latent position t',y='Arc length v',title='E. The coordinate transformation',
 subtitle='Dashed line: constant-speed reference v = L t')+line_theme
# Speeds are compared at corresponding physical locations indexed by original t.
p6 <- ggplot(data.frame(t=t,speed=speed),aes(t,speed))+geom_line(color='#0072B2',linewidth=.85)+
 geom_hline(yintercept=1,color='#D55E00',linetype=2,linewidth=.9)+
 labs(x='Original t (matching locations on the curve)',y='Curve speed',title='F. Speed before and after reparameterization',
 subtitle='Blue: ||df/dt||; orange dashed: ||dg/dv|| = 1')+line_theme
plots <- list(p1,p2,p3,p4,p5,p6)
# Use grid directly to avoid adding a layout dependency to the study library.
draw_all <- function() {
 grid::grid.newpage()
 layout <- grid::grid.layout(4,2,heights=grid::unit(c(.95,.95,.85,.15),'null'))
 grid::pushViewport(grid::viewport(layout=layout))
 for(i in seq_along(plots)) print(plots[[i]],vp=grid::viewport(layout.pos.row=ceiling(i/2),layout.pos.col=(i-1)%%2+1))
 grid::grid.text(paste('Arc length is computed in the selected V13-V14 plane over t in [0,1], using the exact noise-free generating functions.',
 'Feature values and the geometric path are unchanged; only position spacing changes. Hollow circles mark the two endpoints.',
 'If the horizontal coordinate is normalized to w = v/L in [0,1], then ||dg/dw|| = L, not 1.',sep='\n'),
  x=.03,y=.5,just='left',gp=grid::gpar(fontsize=10),vp=grid::viewport(layout.pos.row=4,layout.pos.col=1:2))
 grid::popViewport()
}
png(file.path(out,'true_functions_t_vs_arc_length.png'),width=2400,height=2700,res=180);draw_all();dev.off()
pdf(file.path(out,'true_functions_t_vs_arc_length.pdf'),width=13.333,height=15);draw_all();dev.off()
# A compact first view with matched [0,1] horizontal ranges, explicitly normalized.
p_normalized <- make_function_plot(v_grid/L,g,'Normalized arc length w = v/L',
 'True functions on normalized arc length','Same feature values; constant speed is L on this normalized scale')
ggsave(file.path(out,'normalized_arc_length_functions.png'),p_normalized,width=8,height=4.8,dpi=180,bg='white')
print(checks,row.names=FALSE)
draw_overview <- function() {
 grid::grid.newpage()
 grid::pushViewport(grid::viewport(layout=grid::grid.layout(3,2,heights=grid::unit(c(1,1,.12),'null'))))
 for(i in 1:4)print(plots[[i]],vp=grid::viewport(layout.pos.row=ceiling(i/2),layout.pos.col=(i-1)%%2+1))
 grid::grid.text(paste('The exact same true curve, with a different coordinate: feature values and geometric shape are unchanged.',
 'Arc length uses these two features only. On w = v/L in [0,1], the speed is L, not 1.',sep='\n'),
 x=.03,y=.5,just='left',gp=grid::gpar(fontsize=10),vp=grid::viewport(layout.pos.row=3,layout.pos.col=1:2))
 grid::popViewport()
}
png(file.path(out,'transformation_overview.png'),width=2400,height=1900,res=180);draw_overview();dev.off()
pdf(file.path(out,'transformation_overview.pdf'),width=13.333,height=10.556);draw_overview();dev.off()
