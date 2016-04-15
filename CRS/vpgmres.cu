#include "vpgmres.h"

int VPGMRES_CRS(double *val, int *col, int *ptr, double *bvec, double *xvec, int ndata, int nnz, double eps, int i_max, int rs)
{
  int i,j,k;
  FILE *p_x,*p_his;
  double *rvec, *axvec, *evec, *vvec, *vmtx, *hmtx, *yvec, *wvec, *avvec, *hvvec, *cvec, *svec, *x0vec, *tmpvec, *zmtx, *zvec, *x_0;
  double alpha,wv_ip;
  int count = 0;
  double tmp,tmp2,eps_now=0.0,b_norm;
  int flag = 0;
  int flag_break = 0;
  int error_message, t_error;

  

  rvec = (double *)malloc(sizeof(double)*ndata);
  axvec = (double *)malloc(sizeof(double)*ndata);
  evec = (double *)malloc(sizeof(double)*rs);
  vvec = (double *)malloc(sizeof(double)*ndata);
  vmtx = (double *)malloc(sizeof(double)*ndata*(rs+1));
  hmtx = (double *)malloc(sizeof(double)*ndata*(rs+1));
  yvec = (double *)malloc(sizeof(double)*rs);
  wvec = (double *)malloc(sizeof(double)*ndata);
  avvec = (double *)malloc(sizeof(double)*ndata);
  hvvec = (double *)malloc(sizeof(double)*rs*(rs+1));
  cvec = (double *)malloc(sizeof(double)*rs);
  svec = (double *)malloc(sizeof(double)*rs);
  x0vec = (double *)malloc(sizeof(double)*ndata);
  tmpvec = (double *)malloc(sizeof(double)*ndata);
  zmtx = (double *)malloc(sizeof(double)*ndata*(rs+1));
  zvec = (double *)malloc(sizeof(double)*ndata);
  x_0 = (double *)malloc(sizeof(double)*ndata);



  /* if((p_x = fopen("./output/x.txt", "w")) == NULL){ */
  /*   printf("x.txt open error\n"); */
  /*   exit(1); */
  /* } */
  /* if((p_his = fopen("./output/his.txt", "w")) == NULL){ */
  /*   printf("his.txt open error\n"); */
  /*   exit(1); */
  /* } */

  p_x = FileInit("./output/VPGMRES_x.txt", "w");
  p_his = FileInit("./output/VPGMRES_his.txt", "w");


  for(i=0;i<ndata;i++){
    rvec[i]=0.0;
    axvec[i]=0.0;
    vvec[i]=0.0;
    wvec[i]=0.0;
    avvec[i]=0.0;
    x0vec[i]=0.0;
    xvec[i]=0.0;
    tmpvec[i]=0.0;
    zvec[i]=0.0;
  }

  for(i=0;i<ndata*(rs+1);i++){
    vmtx[i]=0.0;
    hmtx[i]=0.0;
    zmtx[i]=0.0;
  }

  for(i=0;i<rs;i++){
    yvec[i]=0.0;
    cvec[i]=0.0;
    svec[i]=0.0;
    evec[i]=0.0;
  }

  for(i=0;i<rs*(rs+1);i++){
    hvvec[i]=0.0;
  }
  /* b_norm = vector_norm_2(bvec,ndata); */
  b_norm = Double2Norm(bvec, ndata);
  DoubleVecCopy(x_0, xvec, ndata);

  //outer loop 
  for(count=0;count<i_max;){
    //Ax0
    tmp=0.0;
    #pragma omp parallel for private(j) reduction(+:tmp) schedule(static) firstprivate(axvec, val, xvec) lastprivate(axvec)
    for(i=0;i<ndata;i++){
      tmp=0.0;
      for(j=ptr[i];j<ptr[i+1];j++){
        tmp+=val[j] * xvec[col[j]];
      }
      axvec[i]=tmp;
    }
    //r0 = b - Ax0
    for(i=0;i<ndata;i++){
      rvec[i] = bvec[i] - axvec[i];
    }
    //2norm(r)
    /* tmp = vector_norm_2(rvec,ndata); */
    tmp = Double2Norm(rvec, ndata);

    //v0 = r0 / 2norm(r)
    //v0 >> vmtx[0][]
    for(i=0;i<ndata;i++){
      vvec[i] = rvec[i] / tmp;
      vmtx[0*ndata+i] = vvec[i];
    }


    evec[0] = tmp;
    for(i=1;i<rs;i++){
      evec[i]=0.0;
    }

    //inner loop k-> restart
    for(k=0;k<rs-1;k++){
      eps_now = fabs(evec[k]) / b_norm;
      if(verbose)
        printf("%d %.12e\n",count+1, eps_now);
      if(count >= i_max){
        flag_break = 1;
        break;
      }
      fprintf(p_his,"%d %.12e\n",count+1,eps_now);
      //if over eps break
      if(eps_now <= eps){
        solve_Hye(hmtx,yvec,evec,k,ndata);

        //epsilon yv
        for(i=0;i<ndata;i++){
          tmpvec[i]=0.0;
        }
        for(i=0;i<k;i++){
          for(j=0;j<ndata;j++){
            /* tmpvec[j] += yvec[i] * vmtx[i*ndata+j]; */
            tmpvec[j] += yvec[i] * zmtx[i*ndata+j];
          }
        }
        
        //x = x0 + epsilon yv
        for(i=0;i<ndata;i++){
          xvec[i] = x0vec[i] + tmpvec[i];
        }
        flag = 1;
        /* printf("%.12e\n",eps_now); */
        break;
      }
      //TODO
      //inner solver
      /* JOR_CRS(val, col, ptr, vvec, zvec, ndata, inner_eps, inner_i_max, inner_omega); */
      /* cgm_CRS(val, col, ptr, vvec, zvec, ndata, inner_eps, inner_i_max, 0); */
      error_message=InnerSolverSelecter(val, col, ptr, vvec, zvec, ndata, nnz, eps_inner, loop_inner, kskip_inner, fix_inner);
      if(error_message!=0){
        printf("error in inner solver\n");
        return -1;
      }

      /* //Av & W */
      #pragma omp parallel for private(j) reduction(+:tmp) schedule(static) firstprivate(wvec, val, zvec) lastprivate(wvec)
      for(i=0;i<ndata;i++){
        tmp=0.0;
        for(j=ptr[i];j<ptr[i+1];j++){
          /* tmp+=val[j] * vmtx[k*ndata+col[j]]; */
          tmp += val[j] * zvec[col[j]];
        }
        /* avvec[i]=tmp; */
        /* wvec[i]=avvec[i]; */
        wvec[i]=tmp;
      }
      for(i=0;i<ndata;i++){
        zmtx[k*ndata+i] = zvec[i];
      }
      //h_i_k & W  update
      for(i=0;i<=k;i++){
        for(j=0;j<ndata;j++){
          tmpvec[j] = vmtx[i*ndata+j];
        }
      }
      for(i=0;i<=k;i++){
        wv_ip=0.0;
        for(j=0;j<ndata;j++){
          wv_ip+=wvec[j] * vmtx[i*ndata+j];
        }
        hmtx[i*ndata+k] = wv_ip;
        for(j=0;j<ndata;j++){
          wvec[j]=wvec[j] - wv_ip * vmtx[i*ndata+j];
        }
      }
      //h_k+1_k update
      /* tmp=vector_norm_2(wvec,ndata); */
      tmp = Double2Norm(wvec, ndata);
      hmtx[(k+1)*ndata+k]=tmp;

      //v update
      for(i=0;i<ndata;i++){
        vvec[i] = wvec[i] / tmp;
        vmtx[(k+1)*ndata+i] = vvec[i];
      }

      //h_ update
      for(i=0;i<=(k-1);i++){
        tmp=hmtx[i*ndata+k];
        tmp2=hmtx[(i+1)*ndata+k];
        hmtx[i*ndata+k] = cvec[i] * tmp - svec[i] * tmp2;
        hmtx[(i+1)*ndata+k] = svec[i] * tmp + cvec[i] * tmp2;
      }

      //alpha = root(h_kk * h_kk + h_k+1_k * h_k+1_k)
      alpha = sqrt(hmtx[k*ndata+k] * hmtx[k*ndata+k] + hmtx[(k+1)*ndata+k] * hmtx[(k+1)*ndata+k]);

      cvec[k] = hmtx[k*ndata+k] / alpha;
      svec[k] = -hmtx[(k+1)*ndata+k] / alpha;
      evec[k+1] = svec[k] * evec[k];
      evec[k] = cvec[k] * evec[k];
      hmtx[k*ndata+k] = cvec[k] * hmtx[k*ndata+k] - svec[k] * hmtx[(k+1)*ndata+k];
      hmtx[(k+1)*ndata+k] = 0.0;


      count ++;
    }

    if(flag == 1){
      break;
    }
    if(flag_break == 1){
      break;
    }    
    solve_Hye(hmtx, yvec, evec, rs-1, ndata);

    for(i=0;i<ndata;i++){
      tmpvec[i]=0.0;
    }
    for(i=0;i<rs;i++){
      for(j=0;j<ndata;j++){
        /* tmpvec[j] += yvec[i] * vmtx[i*ndata+j]; */
        tmpvec[j] += yvec[i] * zmtx[i*ndata+j];
      }
    }

    for(i=0;i<ndata;i++){
      xvec[i] = x0vec[i] + tmpvec[i];
    }

    for(i=0;i<ndata;i++){
      x0vec[i] = xvec[i];
    }
  }

/*   if(flag == 1){ */
/*     printf("good\n"); */
/*   }else{ */
/*     printf("bad\n"); */
/*   } */
/*   printf("imax = %d roop = %d\n", i_max, count); */
/*   printf("residual = %.12e\n",eps_now); */
/* #ifdef EBUG */
/*   for(i=0;i<ndata;i++){ */
/*     fprintf(p_x, "%d %.12e\n", i, xvec[i]); */
/*   } */
/* #endif */
/* #ifdef TIME */
/*   printf("ElapsedTime = %.6f ms = %.6f s\n", (et-st)*1000, (et-st)); */
/* #endif */
  FileOutPutVec(p_x, xvec, ndata);
  t_error=error_check_CRS(val, col, ptr, bvec, xvec, x_0, ndata);
  printf("|b-ax|2/|b|2=%.1f\n", t_error);

  free(rvec);
  free(axvec);
  free(evec);
  free(vvec);
  free(vmtx);
  free(hmtx);
  free(yvec);
  free(wvec);
  free(avvec);
  free(hvvec);
  free(cvec);
  free(svec);
  free(x0vec);
  free(tmpvec);
  free(zmtx);
  free(zvec);
  free(x_0);

  FileClose(p_x);
  FileClose(p_his);

  if(flag){
    return 1;
  }
  return 2;

}
