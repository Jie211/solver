#include "header.h"

  void
gmres_CRS(double *val,
    const int *col,
    const int *ptr,
    double *bvec,
    double *xvec,
    const int ndata,
    const double eps,
    const int i_max,
    const int rs)
{
  int i,j,k;
  FILE *p_x,*p_his,*p_loop;
  double *rvec, *axvec, *evec, *vvec, *vmtx, *hmtx, *yvec, *wvec, *avvec, *hvvec, *cvec, *svec, *x0vec, *tmpvec;
  double alpha,wv_ip;
  int count = 0;
  double tmp,tmp2,eps_now=0.0,b_norm;
  int flag = 0;
  /* #ifdef TIME */
  /*   double st,et; */
  /* #endif */
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



  if((p_x = fopen("./output/GMRES_x_out.txt", "w")) == NULL){
    printf("x.txt open error\n");
    exit(1);
  }
  if((p_his = fopen("./output/GMRES_his_out.txt", "a")) == NULL){
    printf("his.txt open error\n");
    exit(1);
  }
  if((p_loop = fopen("./output/GMRES_loop_out.txt", "a")) == NULL){
    printf("loop.txt open error\n");
    exit(1);
  }

  for(i=0;i<ndata;i++){
    rvec[i]=0.0;
    axvec[i]=0.0;
    vvec[i]=0.0;
    wvec[i]=0.0;
    avvec[i]=0.0;
    x0vec[i]=0.0;
    xvec[i]=0.0;
    tmpvec[i]=0.0;
  }

  for(i=0;i<ndata*(rs+1);i++){
    vmtx[i]=0.0;
    hmtx[i]=0.0;
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
  /* #ifdef TIME */
  /*   st=gettimeofday_sec(); */
  /* #endif */
  b_norm = vector_norm_2(bvec,ndata);

  //outer loop 
  for(count=0;count<i_max;){
    /* printf("##########on %d##########\n", count); */
    //Ax0
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
    // display_vec(rvec,ndata);
    //2norm(r)
    tmp = vector_norm_2(rvec,ndata);
    // printf("r norm %d %.12e\n", count, tmp);

    //v0 = r0 / 2norm(r)
    //v0 >> vmtx[0][]
    for(i=0;i<ndata;i++){
      vvec[i] = rvec[i] / tmp;
      vmtx[0*ndata+i] = vvec[i];
    }

    // printf("-----------------------\n");
    // for(i=0;i<rs+1;i++){
    //   for(j=0;j<ndata;j++){
    //     printf("%d %d %.12e\n",i ,j ,vmtx[i*ndata+j]);
    //   }
    //   putchar('\n');
    // }
    // printf("-----------------------\n");


    evec[0] = tmp;
    for(i=1;i<rs;i++){
      evec[i]=0.0;
    }

    //inner loop k-> restart
    for(k=0;k<rs-1;k++){
      /* printf("k=%d e_%d=%.12e \n",k,k,evec[k]); */
      eps_now = fabs(evec[k]) / b_norm;
      /* printf("%d %.12e\n", count+1,eps_now); */
#ifdef EBUG
      fprintf(p_his,"%d %.12e\n",count+1,eps_now);
#endif
      //if over eps break
      if(eps_now <= eps){
        solve_Hye(hmtx,yvec,evec,k,ndata);

        //epsilon yv
        for(i=0;i<ndata;i++){
          tmpvec[i]=0.0;
        }
        for(i=0;i<k;i++){
          for(j=0;j<ndata;j++){
            tmpvec[j] += yvec[i] * vmtx[i*ndata+j];
          }
        }
        //x = x0 + epsilon yv
        for(i=0;i<ndata;i++){
          xvec[i] = x0vec[i] + tmpvec[i];
        }
        flag = 1;
        break;
      }

      //Av & W
      for(i=0;i<ndata;i++){
        tmp=0.0;
        for(j=ptr[i];j<ptr[i+1];j++){
          tmp+=val[j] * vmtx[k*ndata+col[j]];
        }
        avvec[i]=tmp;
        wvec[i]=avvec[i];
      }
      /* printf("-----------------------\n"); */
      /* for(i=0;i<rs+1;i++){ */
      /*   for(j=0;j<ndata;j++){ */
      /*     printf("%d %d %.12e\n",i ,j ,vmtx[i*ndata+j]); */
      /*   } */
      /*   putchar('\n'); */
      /* } */
      /* printf("-----------------------\n"); */

      //h_i_k & W  update
      for(i=0;i<=k;i++){
        for(j=0;j<ndata;j++){
          tmpvec[j] = vmtx[i*ndata+j];
        }
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
      tmp=vector_norm_2(wvec,ndata);
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
      /* printf("%.12e %.12e\n", hmtx[k*ndata+k], hmtx[(k+1)*ndata+k]); */
      /* printf("alpha=%.12e\n",alpha); */

      cvec[k] = hmtx[k*ndata+k] / alpha;
      svec[k] = -hmtx[(k+1)*ndata+k] / alpha;
      evec[k+1] = svec[k] * evec[k];
      evec[k] = cvec[k] * evec[k];
      hmtx[k*ndata+k] = cvec[k] * hmtx[k*ndata+k] - svec[k] * hmtx[(k+1)*ndata+k];
      hmtx[(k+1)*ndata+k] = 0.0;


      count ++;
    }

    //printf(" out inner loop\n");
    if(flag == 1){
      break;
    }
    //printf(" cross break check\n");
    solve_Hye(hmtx, yvec, evec, rs-1, ndata);

    for(i=0;i<ndata;i++){
      tmpvec[i]=0.0;
    }
    for(i=0;i<rs;i++){
      for(j=0;j<ndata;j++){
        tmpvec[j] += yvec[i] * vmtx[i*ndata+j];
      }
    }

    for(i=0;i<ndata;i++){
      xvec[i] = x0vec[i] + tmpvec[i];
    }

    // printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    // display_vec(yvec, ndata);
    // printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");

    for(i=0;i<ndata;i++){
      x0vec[i] = xvec[i];
    }
    //printf("-- x update over --\n");
  }
  /* #ifdef TIME */
  /*   et=gettimeofday_sec(); */
  /* #endif */

  /* printf("imax = %d roop = %d\n", i_max, count); */
  /* printf("residual = %.12e\n",eps_now); */
  fprintf(p_his, "\n\n");
  fprintf(p_loop, "%d\n",count+1);
#ifdef EBUG
  for(i=0;i<ndata;i++){
    fprintf(p_x, "%d %.12e\n", i, xvec[i]);
  }
#endif
  /* #ifdef TIME */
  /*   printf("ElapsedTime = %.6f ms = %.6f s\n", (et-st)*1000, (et-st)); */
  /* #endif */
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

  fclose(p_x);
  fclose(p_his);
  fclose(p_loop);
}
