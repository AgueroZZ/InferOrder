source("experiments/estimate_intrinsic_m_smooth_v032/common.R")
paths <- list.files(file.path(study_dir,"results"),pattern="\\.rds$",full.names=TRUE)
rows <- if(length(paths)) do.call(rbind,lapply(paths,function(p)readRDS(p)$row)) else NULL
by_method <- lapply(c("adaptive","forward"),function(method){
  r <- if(!is.null(rows))rows[rows$method==method,,drop=FALSE]else NULL
  list(method=method,completed=if(is.null(r))0L else nrow(r),
       success=if(is.null(r))0L else sum(r$status=="success"),
       exact=if(is.null(r))0L else sum(r$status=="success" & r$estimated_M==r$true_M,na.rm=TRUE),
       under=if(is.null(r))0L else sum(r$status=="success" & r$estimated_M<r$true_M,na.rm=TRUE),
       over=if(is.null(r))0L else sum(r$status=="success" & r$estimated_M>r$true_M,na.rm=TRUE),
       unresolved=if(is.null(r))0L else sum(r$status!="success"))
})
state <- list(updated_at=format(Sys.time(),tz="UTC",format="%Y-%m-%dT%H:%M:%SZ"),
              planned_outcomes=180L,completed_outcomes=length(paths),
              saved_candidates=length(list.files(file.path(study_dir,"candidates"),pattern="\\.rds$")),
              methods=by_method,design_hash=design_hash)
jsonlite::write_json(state,file.path(study_dir,"run_status.json"),pretty=TRUE,auto_unbox=TRUE)
cat(jsonlite::toJSON(state,pretty=TRUE,auto_unbox=TRUE),"\n")
