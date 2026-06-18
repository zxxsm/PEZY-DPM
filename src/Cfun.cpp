/*
 * File: Cfun.cpp
 * Author: Zhuxiaoxong
 * Email: zhuxiaoxiong_cn@163.com
 * Description: 支持多进程MPI/多线程OpenMP/SC3s异构的 DPM
 */
/**********************************************************************/
#ifdef MPI_P
#include <mpi.h>
#endif

#ifdef OMP_P
#include <omp.h>
#endif

#include <vector>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "Cdglbvar.hpp"
#include "Cfun.hpp"
#define DBG_OUT_ON 1

Tcallcnt *callcnt;

/**********************************************************************/
void timers(double *et)
{
	struct timeval t;
	gettimeofday( &t, (struct timezone *)0 );
	*et = t.tv_sec + t.tv_usec*1.0e-6;
}

/**********************************************************************
 * 
 *  Function: called in Fortran.
 *
 **********************************************************************/
extern "C" void sim_(dTdpmsim *dpmsim,	  // s1-/dpmsim/
					 dTdpmsrc *dpmsrc,	  // s2-/dpmsrc/
					 dTdpmmat *dpmmat,	  // s3-/dpmmat/
					 dTdpmscsr *dpmscsr,	  // s4-/dpmscsr/
					 dTdpmscpw *dpmscpw,	  // s5-/dpmscpw/
					 dTxdpmbw *xdpmbw,	  // s6-/xdpmbw/
					 Txdpmq2D *xdpmq2D,	  // s7-/xdpmq2D/
					 dTdpmlam *dpmlam,	  // s8-/dpmlam/
					 dTdpmlbr *dpmlbr,	  // s9-/dpmlbr/
					 dTdpmstpw *dpmstpw,	  // s10-/dpmstpw/
					 dTxbre *xbre,		  // s11-/xbre/
					 dTdpmlph *dpmlph,	  // s12-/dpmlph/
					 dTdpmcmp *dpmcmp,	  // s13-/dpmcmp/
					 dTdpmpap *dpmpap,	  // s15-/dpmpap/
					 dTdpmvox *dpmvox,	  // s16-/dpmvox/
					 Tdpmrpt *dpmrpt,	  // s17-/dpmrpt/
					 dTcutoff *cutoff,	  // s18-/cutoff/
					 dTdpmwck *dpmwck,	  // s19-/dpmwck/
					 dTdpmsub *dpmsub,	  // s20-/dpmsub/
					 Trseed *rseed0,	  // s21-/RSEED/
					 Tsd *sd0,			  // s22-/sd/
					 Tmpi_para *mpi_para, // s23-/mpi_para/
					 Tcallcnt *callcnt0,  // /callcnt/
					 dTdpmesc *dpmesc,	  // /dpmesc/ 主进程中用于传回到FOTRAN打印
					 double *tlaps)		  // lapstime
{
	callcnt=callcnt0;

	//访存优化连续且是单精度
	Trdpmscsr *rfdpmscsr;
	Trdpmscpw *rfdpmscpw;
	Trxdpmbw *rfxdpmbw;
	Trdpmlam *rfdpmlam;
	Trdpmlbr *rfdpmlbr;
	Trdpmstpw *rfdpmstpw;
	Trdpmlph *rfdpmlph;
	Trdpmwck *rfdpmwck;
	Trdpmpap *rfdpmpap;
	Trdpmcmp *rfdpmcmp;
	
	//双精度转换成单精度
	Tdpmsim *fdpmsim;	  // s1-/dpmsim/
	Tdpmsrc *fdpmsrc;	  // s2-/dpmsrc/
	Tdpmmat *fdpmmat;	  // s3-/dpmmat/
	Txbre   *fxbre;		  // s11-/xbre/
	Tdpmvox *fdpmvox;	  // s16-/dpmvox/
	Tcutoff *fcutoff;	  // s18-/cutoff/
	Tdpmsub *fdpmsub;	  // s20-/dpmsub/
	
	//转换成了访存连续且是单精度
	reconstruct(fdpmsim, dpmsim,fdpmsrc, dpmsrc,
				fdpmmat, dpmmat,fxbre, xbre,
				fdpmvox, dpmvox,fcutoff, cutoff,
				fdpmsub, dpmsub,rfdpmscsr, dpmscsr,
				rfdpmscpw, dpmscpw,rfxdpmbw, xdpmbw,
				rfdpmlam, dpmlam,rfdpmlbr, dpmlbr,
				rfdpmstpw, dpmstpw,rfdpmlph, dpmlph,
				rfdpmwck, dpmwck,rfdpmpap, dpmpap,
				rfdpmcmp, dpmcmp);
#ifdef PEZY_SCX
	//在SC3s上进行数据传输和计算
	pzc_simulate(fdpmsim, fdpmvox, fdpmsrc,
				 dpmesc, rfdpmwck, rfdpmcmp,
				 dpmrpt, rfdpmlph, rfdpmpap,
				 fdpmmat, fdpmsub, rfdpmscsr,
				 rfdpmlam, rfdpmstpw, rfdpmscpw,
				 rfdpmlbr, rfxdpmbw, xdpmq2D,
				 fcutoff, fxbre, sd0, rseed0,
				 mpi_para,tlaps);
#else
	//在CPU上进行计算
	simulate(dpmsim, fdpmvox, fdpmsrc,
			 dpmesc, rfdpmwck, rfdpmcmp,
			 dpmrpt, rfdpmlph, rfdpmpap,
			 fdpmmat, fdpmsub, rfdpmscsr,
			 rfdpmlam, rfdpmstpw, rfdpmscpw,
			 rfdpmlbr, rfxdpmbw, xdpmq2D,
			 fcutoff, fxbre,sd0,rseed0,
			 mpi_para,tlaps);	
#endif
}	
/**********************************************************************
 * 
 *  数据进行重构；
 *  
 *
 **********************************************************************/
void reconstruct(Tdpmsim *&fdpmsim, dTdpmsim *dpmsim,
				 Tdpmsrc *&fdpmsrc, dTdpmsrc *dpmsrc,
				 Tdpmmat *&fdpmmat, dTdpmmat *dpmmat,
				 Txbre *&fxbre, dTxbre *xbre,
				 Tdpmvox *&fdpmvox, dTdpmvox *dpmvox,
				 Tcutoff *&fcutoff, dTcutoff *cutoff,
				 Tdpmsub *&fdpmsub, dTdpmsub *dpmsub,
				 Trdpmscsr *&rdpmscsr, dTdpmscsr *dpmscsr,
				 Trdpmscpw *&rdpmscpw, dTdpmscpw *dpmscpw,
				 Trxdpmbw *&rxdpmbw, dTxdpmbw *xdpmbw,
				 Trdpmlam *&rdpmlam, dTdpmlam *dpmlam,
				 Trdpmlbr *&rdpmlbr, dTdpmlbr *dpmlbr,
				 Trdpmstpw *&rdpmstpw, dTdpmstpw *dpmstpw,
				 Trdpmlph *&rdpmlph, dTdpmlph *dpmlph,
				 Trdpmwck *&rdpmwck, dTdpmwck *dpmwck,
				 Trdpmpap *&rdpmpap, dTdpmpap *dpmpap,
				 Trdpmcmp *&rdpmcmp, dTdpmcmp *dpmcmp)
{
	int32 i,j;
	int32 nSize;
	
	/****************************单精度转换************************/
	//s1-/dpmsim/
	nSize=sizeof(Tdpmsim); 
	fdpmsim = (Tdpmsim *) malloc(nSize);
	memset(fdpmsim, 0, nSize);
	fdpmsim->atime = dpmsim->atime;
	fdpmsim->maxhis = dpmsim->maxhis;
	
	//s2-/dpmsrc/
	nSize=sizeof(Tdpmsrc);
	fdpmsrc = (Tdpmsrc *) malloc(nSize);
	memset(fdpmsrc, 0, nSize);
	fdpmsrc->esrc = dpmsrc->esrc;
	fdpmsrc->eabs = dpmsrc->eabs;
	fdpmsrc->eabsph = dpmsrc->eabsph;
	fdpmsrc->param = dpmsrc->param;
	fdpmsrc->intype = dpmsrc->intype;

	//s3-/dpmmat/
	nSize=sizeof(Tdpmmat); //Tdpmmat *Gdpmmat; //s-07
	fdpmmat = (Tdpmmat *) malloc(nSize);
	memset(fdpmmat, 0, nSize);
	for(i=0; i<maxmat; i++)
	{
		fdpmmat->matden[i] = dpmmat->matden[i];
		fdpmmat->zmass[i] = dpmmat->zmass[i];
		fdpmmat->z2mass[i] = dpmmat->z2mass[i];
	}
	fdpmmat->wcc = dpmmat->wcc;
	fdpmmat->wcb = dpmmat->wcb;
	fdpmmat->nmat = dpmmat->nmat;

	//s11-/xbre/
	nSize=sizeof(Txbre); //Txbre *Gxbre; //s-18
	fxbre = (Txbre *) malloc(nSize);
	memset(fxbre, 0, nSize);
	for(i=0; i<maxmat*3; i++)
	{
		fxbre->f0[i] = xbre->f0[i];
	}
	for(i=0; i<maxmat; i++)
	{
		fxbre->bcb[i] = xbre->bcb[i];
	}	

	//s16-/dpmvox/
	nSize=sizeof(Tdpmvox); 
	fdpmvox = (Tdpmvox *) malloc(nSize);
	memset(fdpmvox, 0, nSize);
	fdpmvox->ddx = dpmvox->ddx;
	fdpmvox->ddy = dpmvox->ddy;
	fdpmvox->ddz = dpmvox->ddz;
	fdpmvox->xmid = dpmvox->xmid;
	fdpmvox->ymid = dpmvox->ymid;
	fdpmvox->zmid = dpmvox->zmid;
	for(i=0;i<nxyz;i++)
	{
		fdpmvox->dens[i] = dpmvox->dens[i];
		fdpmvox->mat[i] = dpmvox->mat[i];
	}
	fdpmvox->Unxvox = dpmvox->Unxvox;
	fdpmvox->Unyvox = dpmvox->Unyvox;
	fdpmvox->Unzvox = dpmvox->Unzvox;

	//s18-/cutoff/
	nSize=sizeof(Tcutoff); //Tcutoff *Gcutoff; //s-17
	fcutoff = (Tcutoff *) malloc(nSize);
	memset(fcutoff, 0, nSize);
	fcutoff->wcion = cutoff->wcion;
	fcutoff->wcbre = cutoff->wcbre;	

	//s20-/dpmsub/
	nSize=sizeof(Tdpmsub); //Tdpmsub *Gdpmsub; //s-08
	fdpmsub = (Tdpmsub *) malloc(nSize);
	memset(fdpmsub, 0, nSize);

	fdpmsub->subden = dpmsub->subden;
	fdpmsub->subfac = dpmsub->subfac;
	fdpmsub->substp = dpmsub->substp;
	
	/*********************单精度+数据连续性转换********************/
	// Restruct s4-/dpmscsr/
	nSize = sizeof(Trdpmscsr);
	rdpmscsr = (Trdpmscsr *)malloc(nSize);
	memset(rdpmscsr, 0, nSize);
	for (i = 0; i < nscsr; i++) // Tdpmscsr *dpmscsr; //s-4
	{
		rdpmscsr->scsspx[i * 4] = dpmscsr->scsspa[i];
		rdpmscsr->scsspx[i * 4 + 1] = dpmscsr->scsspb[i];
		rdpmscsr->scsspx[i * 4 + 2] = dpmscsr->scsspc[i];
		rdpmscsr->scsspx[i * 4 + 3] = dpmscsr->scsspd[i];
	}
	rdpmscsr->escsr = dpmscsr->escsr[0];
	rdpmscsr->idless = dpmscsr->idless;
	
	// Restruct s5-/dpmscpw/ 
	nSize = sizeof(Trdpmscpw);
	rdpmscpw = (Trdpmscpw *)malloc(nSize);
	memset(rdpmscpw, 0, nSize);
	for (j = 0; j < nscp; j++) // Tdpmscpw *dpmscpw; //s-5
	{
		for (i = 0; i < maxmat; i++)
		{
			rdpmscpw->scpspx[j * maxmat * 2 + i * 2] = dpmscpw->scpspa[j * maxmat + i];
			rdpmscpw->scpspx[j * maxmat * 2 + i * 2 + 1] = dpmscpw->scpspb[j * maxmat + i];
		}
	}
	rdpmscpw->escp = dpmscpw->escp[0];
	rdpmscpw->idlesc = dpmscpw->idlesc;
	// Restruct s6-/xdpmbw/ 
	nSize = sizeof(Trxdpmbw);
	rxdpmbw = (Trxdpmbw *)malloc(nSize);
	memset(rxdpmbw, 0, nSize);
	for (i = 0; i < nbw; i++) // Txdpmbw *xdpmbw; //s-6
	{
		rxdpmbw->bwspx[i * 4] = xdpmbw->bwspa[i];
		rxdpmbw->bwspx[i * 4 + 1] = xdpmbw->bwspb[i];
		rxdpmbw->bwspx[i * 4 + 2] = xdpmbw->bwspc[i];
		rxdpmbw->bwspx[i * 4 + 3] = xdpmbw->bwspd[i];
	}
	rxdpmbw->ebw = xdpmbw->ebw[0];
	rxdpmbw->idlebw = xdpmbw->idlebw;

	// Restruct s8-/dpmlam/
	nSize = sizeof(Trdpmlam);
	rdpmlam = (Trdpmlam *)malloc(nSize);
	memset(rdpmlam, 0, nSize);
	for (i = 0; i < nlam; i++) // Tdpmlam *dpmlam; //s-8
	{
		rdpmlam->lamspx[i * 4] = dpmlam->lamspa[i];
		rdpmlam->lamspx[i * 4 + 1] = dpmlam->lamspb[i];
		rdpmlam->lamspx[i * 4 + 2] = dpmlam->lamspc[i];
		rdpmlam->lamspx[i * 4 + 3] = dpmlam->lamspd[i];
	}
	rdpmlam->elam = dpmlam->elam[0];
	rdpmlam->idlela = dpmlam->idlela;

	// Restruct s9-/dpmlbr/
	nSize = sizeof(Trdpmlbr);
	rdpmlbr = (Trdpmlbr *)malloc(nSize);
	memset(rdpmlbr, 0, nSize);
	for (j = 0; j < nlra; j++) // Tdpmlbr *dpmlbr; //s-9
	{
		for (i = 0; i < maxmat; i++)
		{
			rdpmlbr->lraspx[j * maxmat * 2 + i * 2] = dpmlbr->lraspa[j * maxmat + i];
			rdpmlbr->lraspx[j * maxmat * 2 + i * 2 + 1] = dpmlbr->lraspb[j * maxmat + i];
		}
	}
	rdpmlbr->elra = dpmlbr->elra[0];
	rdpmlbr->idlelr = dpmlbr->idlelr;

	// Restruct s10-/dpmstpw/
	nSize = sizeof(Trdpmstpw);
	rdpmstpw = (Trdpmstpw *)malloc(nSize);
	memset(rdpmstpw, 0, nSize);
	for (j = 0; j < nst; j++) // Tdpmstpw *dpmstpw; //s-10
	{
		for (i = 0; i < maxmat; i++)
		{
			rdpmstpw->stspx[j * maxmat * 2 + i * 2] = dpmstpw->stspa[j * maxmat + i];
			rdpmstpw->stspx[j * maxmat * 2 + i * 2 + 1] = dpmstpw->stspb[j * maxmat + i];
		}
	}
	rdpmstpw->est = dpmstpw->est[0];
	rdpmstpw->idlest = dpmstpw->idlest;
	
	// Restruct s12-/dpmlph/
	nSize = sizeof(Trdpmlph);
	rdpmlph = (Trdpmlph *)malloc(nSize);
	memset(rdpmlph, 0, nSize);
	for (j = 0; j < nlaph; j++) // Tdpmlph *dpmlph; //s-05
	{
		for (i = 0; i < maxmat; i++)
		{
			rdpmlph->lamphx[j * maxmat * 4 + i * 4] = dpmlph->lampha[j * maxmat + i];
			rdpmlph->lamphx[j * maxmat * 4 + i * 4 + 1] = dpmlph->lamphb[j * maxmat + i];
			rdpmlph->lamphx[j * maxmat * 4 + i * 4 + 2] = dpmlph->lamphc[j * maxmat + i];
			rdpmlph->lamphx[j * maxmat * 4 + i * 4 + 3] = dpmlph->lamphd[j * maxmat + i];
		}
	}
	rdpmlph->elaph = dpmlph->elaph[0];
	rdpmlph->idleph = dpmlph->idleph;
	
	// Restruct s13-/dpmcmp/
	nSize = sizeof(Trdpmcmp);
	rdpmcmp = (Trdpmcmp *)malloc(nSize);
	memset(rdpmcmp, 0, nSize);
	for (j = 0; j < ncmpt; j++) // Tdpmcmp *dpmcmp; //s-04
	{
		for (i = 0; i < maxmat; i++)
		{
			rdpmcmp->comptx[j * maxmat * 4 + i * 4] = dpmcmp->compta[j * maxmat + i];
			rdpmcmp->comptx[j * maxmat * 4 + i * 4 + 1] = dpmcmp->comptb[j * maxmat + i];
			rdpmcmp->comptx[j * maxmat * 4 + i * 4 + 2] = dpmcmp->comptc[j * maxmat + i];
			rdpmcmp->comptx[j * maxmat * 4 + i * 4 + 3] = dpmcmp->comptd[j * maxmat + i];
		}
	}
	rdpmcmp->ecmpt = dpmcmp->ecmpt[0];
	rdpmcmp->idlecp = dpmcmp->idlecp;

	// Restruct s15-/dpmpap/
	nSize = sizeof(Trdpmpap);
	rdpmpap = (Trdpmpap *)malloc(nSize);
	memset(rdpmpap, 0, nSize);
	for (j = 0; j < npair; j++) // Tdpmpap *dpmpap; //s-15
	{
		for (i = 0; i < maxmat; i++)
		{
			rdpmpap->pairpx[j * maxmat * 4 + i * 4] = dpmpap->pairpa[j * maxmat + i];
			rdpmpap->pairpx[j * maxmat * 4 + i * 4 + 1] = dpmpap->pairpb[j * maxmat + i];
			rdpmpap->pairpx[j * maxmat * 4 + i * 4 + 2] = dpmpap->pairpc[j * maxmat + i];
			rdpmpap->pairpx[j * maxmat * 4 + i * 4 + 3] = dpmpap->pairpd[j * maxmat + i];
		}
	}
	rdpmpap->epair = dpmpap->epair[0];
	rdpmpap->idlepp = dpmpap->idlepp;

	// Restruct s19-/dpmwck/
	nSize=sizeof(Trdpmwck); 
	rdpmwck = (Trdpmwck *) malloc(nSize);
	memset(rdpmwck, 0, nSize);
	for(i=0; i<nwck; i++) // Tdpmwck *dpmwck; //s-19
	{
		rdpmwck->awck[i*2] = dpmwck->a0wck[i];
		rdpmwck->awck[i*2+1] = dpmwck->a1wck[i];
	}
	rdpmwck->idlewk = dpmwck->idlewk;
	rdpmwck->wcke0 = dpmwck->wcke0;
}
/**********************************************************************
 * 
 *  Function: main  simulation.
 *
 **********************************************************************/
void simulate(dTdpmsim *dpmsim, Tdpmvox *dpmvox, Tdpmsrc *dpmsrc,
			  dTdpmesc *dpmesc_pro, Trdpmwck *dpmwck, Trdpmcmp *dpmcmp,
			  Tdpmrpt *dpmrpt, Trdpmlph *dpmlph, Trdpmpap *dpmpap,
			  Tdpmmat *dpmmat, Tdpmsub *dpmsub, Trdpmscsr *dpmscsr,
			  Trdpmlam *dpmlam, Trdpmstpw *dpmstpw, Trdpmscpw *rdpmscpw,
			  Trdpmlbr *dpmlbr, Trxdpmbw *xdpmbw, Txdpmq2D *xdpmq2D,
			  Tcutoff *cutoff, Txbre *xbre, Tsd *sd0,Trseed *rseed0,
			  Tmpi_para *mpi_para, double *tlaps)
{	
	int32 CCNT_SUM_PR=0;//记录总SEED调用次数
	int32 CCNT_SUM_TH=0;//记录多线程调用次数
	
	double time0, time1;

#ifdef OMP_P
	//在使用openmp时，Trdpmstck和Trdpmesc两个变量分配在堆上；与SCx保持一致
	//为了保证访存的不跨越性，Trdpmstck和Trdpmesc都需要扩展至与cacheline对齐
	//dpmstck和dpmesc需要进行线程数的扩展
	int32 err;
	Trdpmstck *Gdpmstck;
	err = posix_memalign((void**)&(Gdpmstck), 0x40, THREAD_NUM*sizeof(Trdpmstck));

	Trdpmesc  *Gdpmesc;
	err = posix_memalign((void**)&(Gdpmesc ), 0x40, THREAD_NUM*sizeof(Trdpmesc));
	memset(Gdpmesc, 0, THREAD_NUM*sizeof(Trdpmesc));

#else
	//dpmstck和dpmesc参与计算
	Trdpmstck *dpmstck=(Trdpmstck *)malloc(sizeof(Trdpmstck));

	//与4有区别，Trdpmesc是新的数据结构
	Trdpmesc *dpmesc = (Trdpmesc*)malloc(sizeof(Trdpmesc));
	memset(dpmesc, 0, sizeof(Trdpmesc));

#endif
#ifdef MPI_P
    MPI_Barrier(MPI_COMM_WORLD);
#endif
	//开始计时
	timers(&time0);

#ifdef OMP_P
//启动openmp
#pragma omp parallel num_threads(THREAD_NUM) reduction(+:CCNT_SUM_TH)
	{
		//获取线程数
		int32 thread_id=omp_get_thread_num();
		//获取全局id
		int32 gloabl_id=mpi_para->myrank*THREAD_NUM+omp_get_thread_num();
		//获取全局核数（进程×线程）
		int32 gloabl_num=mpi_para->pro_num*THREAD_NUM;
		//dpmstck和dpmesc参与计算
		Trdpmstck *dpmstck=&Gdpmstck[thread_id];
		Trdpmesc *dpmesc =&Gdpmesc [thread_id];

#else
		//线程数设为1
		int32 thread_id=0;
		//获取全局id
		int32 gloabl_id=mpi_para->myrank;
		//获取全局核数（进程×线程）
		int32 gloabl_num=mpi_para->pro_num;
#endif
#ifdef DEBUG
#ifdef OMP_P
		if (gloabl_id==0)
		{
			//目标区域的大小
			printf("All values are nxini=%d, nxfin=%d, nyini=%d, nyini=%d, nzini=%d, nzfin=%d\n", nxini, nxfin, nyini, nyfin, nzini, nzfin);
			printf("PRO-NUM:%d,TH-NUM:%d,GL-NUM:%d,PART-NUM-ALL:%d\n",mpi_para->pro_num,THREAD_NUM,gloabl_num,_maxhis);

		}
#else
		if (gloabl_id == 0)
		{
			printf("PRO-NUM:%d,TH-NUM:%d,GL-NUM: 1,PART-NUM-ALL:%d\n", mpi_para->pro_num,gloabl_num, _maxhis);
		}
#endif
#endif
		/*******************计算每个ID的粒子数*********************/	
		int32 nprtcl= Ccomp_np(_maxhis,gloabl_id,gloabl_num);

		//seed 初始化
		Trseed rseed;
		rseed.IICNT = 0;
		if (gloabl_num == 1)//只在串行计算的时候才使用原始DPM的SEED
		{
			rseed.IISEED1 = rseed0->IISEED1; 
			rseed.IISEED2 = rseed0->IISEED2;
		}
		else//其他情况都使用自己生成的SEED
		{
			rseed.IISEED1 = sd0->Iseed[gloabl_id * 2];
			rseed.IISEED2 = sd0->Iseed[gloabl_id * 2 + 1];
		}
#ifdef DEBUG
		printf("PRO-ID:%d,TH-ID:%d,GL-ID:%d,PART-NUM:%d,RAND-NUM:%d %d\n",
		mpi_para->myrank,thread_id,gloabl_id,nprtcl,rseed.IISEED1,rseed.IISEED2);
#ifdef OMP_P
		#pragma omp barrier
#endif
#endif

		/****************定义线程私有变量，并进行初始化，这部分分配在栈里*************/
		// 新的计数器
		Tdpmesc_new dpmesc_new;
		
		// Trdpmstck-nsec=0 对应原始序的init
		_nsec = 0;
		// p1-原始序中在source中直接初始化，这里不要任何操作
		Tdpmpart dpmpart;
		// p2-原始程序中在inigeo中进行初始化,这里一一对应
		Tcgeom3 cgeom3;
		cgeom3.dxx = _ddx;
		cgeom3.dyy = _ddy;
		cgeom3.dzz = _ddz;
		cgeom3.idxx = 1.0 / cgeom3.dxx;
		cgeom3.idyy = 1.0 / cgeom3.dyy;
		cgeom3.idzz = 1.0 / cgeom3.dzz;
		// p3-记录电子发生作用的信息,electr中初始化和计算
		Tdpmjmp dpmjmp;
		// p4-记录电子发生作用的时移动穿越哪个截面
		Tcgo cgo;

		/******************************定义线程结束**********************************/

		/*****************************开始粒子主循环*********************************/
		for (int32 i = 1; i <= nprtcl; i++)
		{
			dpmpart.nnhist = i;
			
			//每个粒子都需要对dpmesc_new初始化为0;
			memset(&dpmesc_new, 0, sizeof(Tdpmesc_new));

			Csource(&dpmpart, dpmvox, dpmsrc, &cgeom3, &rseed);
			do
			{
				if (dpmpart.pptype == 0)
				{
					Cphoton(&dpmpart, dpmstck, dpmvox, dpmsrc, &dpmesc_new, dpmwck,
							dpmcmp, dpmrpt, &cgeom3, dpmlph, dpmpap, &rseed);
				}
				else
				{
					Celectr(&dpmpart, dpmvox, dpmmat, &dpmjmp, dpmsrc, dpmsub, dpmsim,
							dpmscsr, dpmlam, &cgo, &cgeom3, dpmstpw, rdpmscpw, dpmlbr,
							&dpmesc_new, dpmrpt, xdpmbw, xdpmq2D, dpmstck, cutoff, xbre, &rseed);
					if (dpmpart.pptype == 1)
					{
						Cputann(&dpmpart, dpmstck, dpmsrc, &rseed);
					}
				}
			} while (Cscndry(&dpmpart, dpmstck) == 1);
			
			//每个粒子完成计算后后将dpmesc_new——>dpmesc中
			Cdumpe_new(dpmesc,&dpmesc_new);
		}
		
		/*****************************结束粒子主循环*********************************/
		//CCNT规约
		CCNT_SUM_TH+=rseed.IICNT;
#ifdef DEBUG
		printf("PRO-ID:%d,TH-ID:%d,GL-ID:%d,CCNT:%d\n",
		mpi_para->myrank,thread_id,gloabl_id,rseed.IICNT);
#ifdef OMP_P
		#pragma omp barrier
#endif
#endif		
#ifdef OMP_P
	}
#endif

#ifdef OMP_P
	//OMP线程级规约dpmesc
	for (int32 i = 1; i < THREAD_NUM; i++)
	{
#pragma omp parallel for num_threads(THREAD_NUM)
		for (int32 j = 0; j < nxyz; j++)
		{
			Gdpmesc[0].escorex[j*2]   = Gdpmesc[0].escorex[j*2] + Gdpmesc[i].escorex[j*2];
			Gdpmesc[0].escorex[j*2+1] = Gdpmesc[0].escorex[j*2+1] + Gdpmesc[i].escorex[j*2+1];
		}
	}

	//指针dpmesc指向Gdpmesc
	Trdpmesc  *dpmesc =&Gdpmesc[0];

#endif

#ifdef MPI_P
	//MPI进程级规约dpmesc
	Trdpmesc *dpmesc_reduce = (Trdpmesc *)malloc(sizeof(Trdpmesc));
	MPI_Reduce(dpmesc, dpmesc_reduce, 2 * nxyz, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&CCNT_SUM_TH, &CCNT_SUM_PR, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
	CCNT_SUM_PR=CCNT_SUM_TH;
#endif

#ifdef DEBUG	
	if (mpi_para->myrank==0)
		printf("CCNT_SUM:%d\n",CCNT_SUM_PR);
#endif

	timers(&time1);
	
#ifdef DEBUG2
	//验证：每个进程的计算结果是否正确
	if (mpi_para->myrank==0)
	{
		for(int32 i=1; i<=nxyz; i++)
			if (_escore(i)!=0.0)
				printf("i=%d,esc=%lf,esc2=%lf\n",i,_escore(i),_escor2(i));
	}
	// printf("xxxxx,esc=%lf,esc2=%lf\n",_escore(279001),_escor2(279001));
#endif

// 数据还原顺序
#ifdef MPI_P
	//不能操作dpmesc_pro,不然会报错gcc bug.
	//icc不会有问题，可能与FORTRAN 与 C++互相调用有关
	for(int i=0; i<nxyz; i++)
	{
		dpmesc->escorex[i]=dpmesc_reduce->escorex[2*i];
	}
	for(int i=0; i<nxyz; i++)
	{
		dpmesc->escorex[i+nxyz]=dpmesc_reduce->escorex[2*i+1];
	}

	memcpy(dpmesc_pro,dpmesc,sizeof(Trdpmesc));//数据拷贝回到 dpmesc_pro

#else
	for(int32 i=1; i<=nxyz; i++)
	{
		dpmesc_pro->escore[i-1]=_escore(i);
		dpmesc_pro->escor2[i-1]=_escor2(i);
	}
#endif

#ifdef MPI_P
	free(dpmesc_reduce);
#endif

#ifdef OMP_P
	// 释放OMP开辟的内存
	free(Gdpmesc);
	free(Gdpmstck);
#else
	// 释放Trdpmesc
	free(dpmesc);
	free(dpmstck);
#endif

	*tlaps = time1 - time0;
}
//OK
/**********************************************************************
 * 
 *  Function: this is an adapted version of subroutine ranecu written 
 *		by f. james (comput. phys. commun. 60 (1990) 329-344), 
 *		which has been modified to give a single random number at each
 *		call.
 *		the 'seeds' iseed1 and iseed2 must be initialized in the
 *		main program and transferred through the named common block
 *		/rseed/.
 *
 **********************************************************************/
UREAL Crng(Trseed *rseed)
{
	int32 I1, I2, IZ;
	UREAL rng;
	
	ICNT = ICNT + 1;
	I1=ISEED1/53668;
	ISEED1=40014*(ISEED1-I1*53668)-I1*12211;
	if(ISEED1<0) ISEED1=ISEED1+2147483563;

	I2=ISEED2/52774;
	ISEED2=40692*(ISEED2-I2*52774)-I2*3791;
	if(ISEED2<0) ISEED2=ISEED2+2147483399;

	IZ=ISEED1-ISEED2;
	if(IZ<1) IZ=IZ+2147483562;
	rng = IZ*USCALE;
	
	return rng;
}//OK 

/**********************************************************************
 * 
 *  Function: Creates a new primary particle state 
 *
 **********************************************************************/
void Csource(Tdpmpart* dpmpart, Tdpmvox *dpmvox, Tdpmsrc *dpmsrc, 
	Tcgeom3 *cgeom3,Trseed *rseed) //Trseed *rseed;
{
	ptype = _intype;
	energy = _esrc;
	vx = 0.0E0f;
	vy = 0.0E0f;
	vz = 1.0E0f;
	x = _xmid+_param*(Crng(rseed)-0.5E0f);
	y = _ymid+_param*(Crng(rseed)-0.5E0f);
	z = +zero;	//z = +zero
	Cwhere(dpmpart, cgeom3, dpmvox);
}//OK-1

/**********************************************************************
 * 
 *  Function: Locates the particle and sets vox in 3D 
 *	Input:	particle position
 *	Output: absvox -> absolute address of the current voxel
 *
 **********************************************************************/
void Cwhere(Tdpmpart *dpmpart, Tcgeom3 *cgeom3, Tdpmvox *dpmvox)
{
/* 	int32 a, b, c;
	int32 tmp;
	
	a = z * idz + 1;
	b = y * idy + 1;
	c = x * idx + 1;
	
#if DBG_OUT_ON == 0		
	cnt1=cnt1+1;
#endif
	zvox = a;
	if((a<1)||(a>Unzvox))
	{		
		absvox = 0;
	} 
	else if((b<1)||(b>Unyvox))
	{
#if DBG_OUT_ON == 0			
		cnt2=cnt2+1;
#endif
		yvox = b;
		absvox = 0;
	} 
	else if((c<1)||(c>Unxvox))	
	{
#if DBG_OUT_ON == 0	
		cnt2=cnt2+1;
#endif
		yvox = b;
		xvox = c;
		absvox = 0;
	} 
	else
	{
#if DBG_OUT_ON == 0			
		cnt2=cnt2+1;
#endif
		yvox = b;
		xvox = c;
		absvox = a + (b-1)*Unzvox+(c-1)*Unzvox*Unyvox;
#if DBG_OUT_ON == 0			
		cnt3=cnt3+1;
#endif
	}
#if DBG_OUT_ON == 0			
	printf("zvox, yvox, xvox are %12d, %12d, %12d \n", zvox, yvox, xvox);
#endif */
/* 	zvox = z * idz + 1;
	if((zvox<1)||(zvox>Unzvox))
	{		
		absvox = 0;
	} else
	{
		yvox = y*idy+1;
		if((yvox<1)||(yvox>Unyvox))
		{
			absvox = 0;
		} 
		else 
		{
			xvox = x*idx+1;
			if((xvox<1)||(xvox>Unxvox))
			{
				absvox = 0;
			} else
			{
				absvox = zvox+(yvox-1)*Unzvox+(xvox-1)*Unzvox*Unyvox;
			}
		} 	
	} */	

	zvox = z * idz + 1;
	if((zvox<1)||(zvox>_Unzvox))
	{		
		absvox = 0;
		return;
	}
	yvox = y*idy+1;
	if((yvox<1)||(yvox>_Unyvox))
	{
		absvox = 0;
		return;
	}
	xvox = x*idx+1;
	if((xvox<1)||(xvox>_Unxvox))
	{
		absvox = 0;
		return;
	}
	absvox = zvox+(yvox-1)*_Unzvox+(xvox-1)*_Unzvox*_Unyvox;
	return;
}//OK-1 

/**********************************************************************
 * 
 *  Function: Stores the inverse of vz to save time while in cango().
 *		It should be called every time {vx,vy,vz} changes.  
 *
 **********************************************************************/
void Csetv(Tdpmpart* dpmpart, Tcgeom3 *cgeom3)
{
	if(vz != 0.0E0f)
	{
		ivz =  1.0E0f/vz;
	}
	else
	{
		ivz = inf;
	}

	if(vy != 0.0E0f)
	{
		ivy =  1.0E0f/vy;
	}
	else
	{
		ivy = inf;
	}

	if(vx != 0.0E0f)
	{
		ivx =  1.0E0f/vx;
	}
	else
	{
		ivx = inf;
	}
} //OK-2

/**********************************************************************
 * 
 *  Function: Informs about interface crossings in a 3D voxel 
 *		geometry and prepares information to move to next voxel 
 *      计算三维体素几何体中的界面交叉点，（电子在某个体素中，根据各个
 * 	    方向体素面的距离和速度计算出穿越的截面是哪个）
 *      并为移动到下一个体素准备信息 
 *	Input:	
 *      {x,y,z} -> position vector (cm)                            
 *      {vx,vy,vz} -> normalized direction vector                  
 *      vox -> voxel where the particle starts  
 *    Output:                                                      
 *      -> distance to intersection with nearest voxel wall (cm)
 *      与最近体素 面相交的距离（厘米）   
 *      index -> {1,2,3} depending on which walls are intersected 
 * 		决定与哪些体素面相交 
 *      dvox -> {+1,-1} depending on whether part. goes forw-backw
 *      决定是向前或向后移动   
 **********************************************************************/
UREAL Cinters(Tdpmpart *dpmpart, Tcgeom3 *cgeom3, Tcgo *cgo)
{
	UREAL inters, smaybe;
	
//om: Checking out all the voxel walls for the smallest distance...
	if(ivz>0.0E0f)
	{
		inters = (zvox * dz - z)*ivz;
        index = 3;
        dvox = +1;	
	}else
	{
		inters = ((zvox-1) * dz - z)*ivz;
        index = 3;
        dvox = -1;		
	}

	if(ivy>0.0E0f)
	{
		smaybe = (yvox * dy - y)*ivy;
		if(smaybe < inters)
		{
			inters = smaybe;
			index = 2;
			dvox = +1;
		}
	} else
	{
		smaybe = ((yvox-1) * dy - y)*ivy;
		if(smaybe < inters)
		{
			inters = smaybe;
			index = 2;
			dvox = -1;
		}
	}

	if(ivx>0.0E0f)
	{
		smaybe = (xvox * dx - x)*ivx;
		if(smaybe < inters)
		{
			inters = smaybe;
			index = 1;
			dvox = +1;	
		}
	} else
	{
		smaybe = ((xvox-1) * dx - x)*ivx;
		if(smaybe < inters)
		{
			inters = smaybe;
			index = 1;
			dvox = -1;
		}
	}
	
//om: Make sure we won't get neg value to avoid interpretation problems...
	if (inters < 0.0E0f) inters = 0.0E0f;
	return inters;
} //OK-2

/**********************************************************************
 * 
 *  Function: Restricted stopping power --linear interpolation 
 *	Input:	
 *      matid -> material id#                            
 *      e -> energy in eV  --kinetic energy--                   
 *    Output:                                                      
 *      StopPow in eV*cm^2/g
 * 
 **********************************************************************/
UREAL Crstpip(int32 Lmatid, UREAL e, Trdpmstpw *dpmstpw)
{
	int32 i;
	UREAL rstpip;
	
	i = _idlest * (e-_est) + 1;
	rstpip = _stspa(Lmatid,i) + e*_stspb(Lmatid,i);
	
	return rstpip;
}//Ok

/**********************************************************************
 * 
 *  Function: Inverse 1st transport MFP --linear interpolation 
 *	Input:	
 *      matid -> material id#                            
 *      e -> energy in eV  --kinetic energy--                   
 *    Output:                                                      
 *      lambda_1^{-1} in cm^2/g
 * 
 **********************************************************************/
UREAL Cscpwip(int32 Lmatid, UREAL e, Trdpmscpw *dpmscpw)
{
	int32 i;
	UREAL scpwip;
	
	i = _idlesc*(e-_escp)+1;
#if DBG_OUT_ON == 0	
	printf("Thread %d: Lmatid=%d, i=%d, e=%6.4E!\n", 0, Lmatid, i, e);
#endif
	scpwip = _scpspa(Lmatid,i)+e*_scpspb(Lmatid,i);
	
	return scpwip;
} //OK-2

/**********************************************************************
 * 
 *  Function: Bremsstrahlung inverse MFP --linear interpolation 
 *	Input:	
 *      matid -> material id number                            
 *      ie -> -1/energy in eV^-1  --kinetic energy--                   
 *    Output:                                                      
 *      inverse mean free path in cm^2/g 
 * 
 **********************************************************************/
UREAL Cilabip(int32 Lmatid, UREAL ie, Trdpmlbr *dpmlbr)
{
	int32 i;
	UREAL de,ilabip;
	
	de = ie-_elra;
	if(de > 0.0E0f)
	{
        i = _idlelr*de+1;
        ilabip = _lraspa(Lmatid,i)+ie*_lraspb(Lmatid,i);	
	} 
	else
	{
		ilabip = zero;
	}

	return ilabip;
} //OK-2

/**********************************************************************
 * 
 *  Function: Deposites energy in the corresponding counters 
 *	Input:	
 *      edep -> energy being deposited (eV)
 *  Worthing Vectorizing ???
 **********************************************************************/
void Cscore(UREAL edep, Tdpmpart *dpmpart, Tdpmesc_new *dpmesc_new, 
	Tdpmrpt *dpmrpt)
{
//	cnt1 = cnt1 +1;
//     *** Do not score if outside the RoI:
#if DBG_OUT_ON == 0			
		printf("%12d %12d %12d\n", zvox, yvox, xvox);
#endif
	// printf("lasthi(absvox)=%d,nhist=%d\n",lasthi(absvox),nhist);
	if((xvox < nxini) || (xvox > nxfin) || (yvox < nyini) ||
		(yvox > nyfin) || (zvox < nzini) || (zvox > nzfin))
	{
		// nothing
	} else
	{
		for (int32 i=0;i<tail;i++)
		{
			//如果absvox有索引相同的则相加,并且返回
			if (absvox==esidx(i))
			{
				estmp(i)+=edep;
				return;
			}
		}
		estmp(tail)=edep;   //记录值
		esidx(tail)=absvox; //记录位置
		tail++;             //尾部加1
				
		
		if (tail >= maxns)
		{
			printf("Cscore:error: Stack is full, enlarge score maxns");
			exit(1);
		}
	}
} //OK-1

/**********************************************************************
 * 
 *  Function: Transport of e- below the nominal Eabs until absorption 
 *	Input:	
 *      {x,y,z} -> take off location
 *		{vx,vy,vz} -> direction of flight
 *		energy -> initial kinetic energy
 *    Output:                                                      *
 *      escore -> energy deposited by CSDA  
 *
 **********************************************************************/
void Csubabs(Tdpmpart *dpmpart, Tdpmvox *dpmvox, Tdpmsub *dpmsub,
	Tdpmesc_new *dpmesc_new, Tcgo *cgo, Tdpmrpt *dpmrpt, Tcgeom3 *cgeom3)
{
	UREAL de, ivoxden, s;
//       Determine if further transport is needed
	while( 1>0 )
	{
		ivoxden = _dens(absvox);
		if(ivoxden > energy*_subfac) 
		{
			Cscore(energy, dpmpart, dpmesc_new, dpmrpt);
			break;
		}
		s = Cinters(dpmpart, cgeom3, cgo);
		de = _substp*ivoxden*s;
		energy = energy-de;
		Cscore(de, dpmpart, dpmesc_new, dpmrpt);
		Cchvox(dpmpart, cgo, dpmvox);
		if(absvox==0)
		{
			break;
		}
        x = x+s*vx;
        y = y+s*vy;
        z = z+s*vz;
	} 
}//OK-2

/**********************************************************************
 * 
 *  Function: Changes voxel according to the information 
 *		passed by inters() 
 *	Input:	
 *      vox -> voxel where the particle starts
 *		index -> {1,2,3} depending on which walls are intersected 
 *		dvox -> {+1,-1} depending on whether part. goes forw-backw
 *    Output:                                                      
 *      vox -> voxel where the final position lies  
 *
 **********************************************************************/
void Cchvox(Tdpmpart *dpmpart, Tcgo *cgo, Tdpmvox *dpmvox)
{
	if(index ==3 )
	{
		zvox = zvox+dvox;
		if((zvox>0)&&(zvox<=_Unzvox))
		{
			absvox = absvox+dvox;
		} else
		{
			absvox = 0;
		}
	}
	else if(index ==2 )
	{
		yvox = yvox+dvox;
		if((yvox>0)&&(yvox<=_Unyvox))
		{
			absvox = absvox+dvox*_Unzvox;
		} else
		{
			absvox = 0;
		}		
	}
	else
	{
		xvox = xvox+dvox;
		if((xvox>0)&&(xvox<=_Unxvox))
		{
			absvox = absvox+dvox*_Unzvox*_Unyvox;
		} else
		{
			absvox = 0;
		}
	}
} //OK-2

/**********************************************************************
 * 
 *  Function: Transports the particle following a rectiliniar 
 *		trajectory, taking care of interface crossings and keeping 
 * 		track of energy losses in the corresponding counters 
 *		(using CSDA).  
 *	Input:	
 *      e- initial state 
 *      Fuel variables affecting the flight                        
 *      event -> event that stopped flight last time; in addition  
 *                to the output codes,                             
 *                -1 new particle                                  
 *                20 end of elastic step  
 *    Output:                                                      
 *      e- final state                                             
 *      Remaining fuel vars                                        
 *      event -> kind of event that causes flight to stop:         
 *                1 run out of energy; absorbed                    
 *                2 run out of elastic fuel                        
 *                3 run out of Moller fuel                         
 *                4 run out of bremsstrahlung fuel                 
 *               99 escaped from universe 
 * 
 ***********************************************************************/
void Cflight(Tdpmpart *dpmpart, Tdpmvox *dpmvox, Tdpmmat *dpmmat, 
	Tdpmjmp *dpmjmp, Tdpmsrc *dpmsrc, Tdpmsub *dpmsub, Tcgo *cgo,
	Tcgeom3 *cgeom3, Trdpmstpw *dpmstpw, Trdpmscpw *dpmscpw, 
	Trdpmlbr *dpmlbr, Tdpmesc_new *dpmesc_new, Tdpmrpt *dpmrpt)
{
	UREAL sback, s, dedx, newe, news;
	UREAL de, infuel;

//     *** Init according to the dynamic vars changed:
	if((event==2) || (event==-1))
	{
		Csetv(dpmpart, cgeom3);
		smax = Cinters(dpmpart, cgeom3, cgo);
	}
	event = 0;
	
//     *** Loop until it runs out of fuel:	
	while(1>0){
//		printf("It's flighting, event: %12d\n", event);
//		printf("ISEED1,ISEED2,ICNT: %12d %12d %12d\n", ISEED1,ISEED2,ICNT);
		s = smax;

//       *** Calculate fuel burn rate in the current voxel:	
		matid = _mat(absvox);
		voxden = _dens(absvox);
		dedx = Crstpip(matid,energy,dpmstpw)*voxden;
		newe = Cmax((energy-0.5f*dedx*s),_eabs);
		burnel = Cscpwip(matid,newe,dpmscpw)*voxden;
		burnmo = _zmass(matid)*voxden;
		burnbr = voxden*Cilabip(matid,-1.0E0f/newe,dpmlbr);

//		*** Burn Moller fuel:
        fuelmo = fuelmo-s*burnmo;
        if (fuelmo<0.0E0f)
		{
			sback = -fuelmo/burnmo;
			s = s-sback;
			fuelmo = 0.0E0f;
			event = 3;
        }
		
//       *** Burn bremss fuel:
		fuelbr = fuelbr-s*burnbr;
		if (fuelbr<0.0E0f)
		{
			sback = -fuelbr/burnbr;
			s = s-sback;
			fuelmo = fuelmo+sback*burnmo;
			fuelbr = 0.0E0f;
			event = 4;
		}
		
//       *** Burn elastic fuel:
        infuel = fuelel;
        fuelel = fuelel-s*burnel;
		if (fuelel<0.0E0f)
		{
//         *** Refine calculation of scattering 1st MFP:
			news = infuel/(Cscpwip(matid,energy,dpmscpw)*voxden);
			newe = Cmax((energy-0.5f*dedx*news),_eabs);
			news = infuel/(Cscpwip(matid,newe,dpmscpw)*voxden);
			if (news>s) news = s;
			sback = s-news;
			s = news;
			fuelmo = fuelmo+sback*burnmo;
			fuelbr = fuelbr+sback*burnbr;
			fuelel = 0.0E0f;
			event = 2;	
		}
		
//       *** Accounting for continous energy loss:
        newe = Cmax((energy-0.5E0f*dedx*s),_eabs);
        de = s*Crstpip(matid,newe,dpmstpw)*voxden;
        energy = energy-de;
        Cscore(de,dpmpart,dpmesc_new,dpmrpt);
		if (energy<_eabs)
		{
//         *** Determine if subEabs transport is needed:
			if(voxden>_subden)
			{
				Cscore(energy,dpmpart,dpmesc_new,dpmrpt);
			} else
			{
				Csubabs(dpmpart, dpmvox, dpmsub, dpmesc_new, cgo, 
					dpmrpt,cgeom3);
			}
			event = 1;
			break;
		}

//       *** Move the electron:
		x = x+s*vx;
		y = y+s*vy;
		z = z+s*vz;
		smax = smax-s;

//       *** Check whether an interaction has not ocurred:
		if (event==0)
		{
			Cchvox(dpmpart, cgo, dpmvox);
			if(absvox==0)
			{
				event = 99;
				break;
			}
			smax = Cinters(dpmpart, cgeom3, cgo);
			continue;
		}
		
//       *** Otherwise, run out of fuel so return:
		break;
	}
} //OK-2

/**********************************************************************
 * 
 *  Function: 3spline interpolation for scattering strength as a	
 *		function of kinetic energy; this quantity is related to 
 *		the step length    
 *	Input:	
 *      e -> kinetic energy in eV 
 *  Output:                                                      
 *      K = scattering strength = integ{ds/lambda1(s)} k
 * 
 ***********************************************************************/
UREAL Cstepip(UREAL e, Trdpmscsr *dpmscsr)
{
	int32 i;
	UREAL stepip;

	i = _idless*(e-_escsr)+1;
	stepip = _scsspa(i)+e*(_scsspb(i)+e*(_scsspc(i)+e*_scsspd(i)));
	
	return stepip;
} //OK-2

/**********************************************************************
 * 
 *  Function: Moller mean free path, 3-spline interpolation     
 *	Input:	
 *      ie -> -1/energy in eV^-1  --kinetic energy-- 
 *  Output:                                                      
 *      mean free path in g/cm^2
 * 
 ***********************************************************************/
UREAL Clamoip(UREAL ie, Trdpmlam *dpmlam)
{
	int32 i;
	UREAL de, lamoip;
	
	de = ie-_elam;
	if (de>0.0E0f)
	{
        i = _idlela*de+1;
        lamoip = 1.0E0f/(_lamspa(i)+ie*(_lamspb(i)+ie*(_lamspc(i)
			+ie*_lamspd(i))));
	}
	else
	{
		lamoip = +inf;
	}
	
	return lamoip;
} //OK-2

/**********************************************************************
 * 
 *  Function: 3spline interpolation for bw as a function energy      
 *	Input:	
 *      ie -> -1/energy in eV^-1  --kinetic energy--
 *  Output:                                                      
 *      bw, broad screening parameter that gets flattest q surf
 *
 ***********************************************************************/
UREAL Cxbwip(UREAL ie, Trxdpmbw *xdpmbw)
{
	UREAL xbwip;
	int32 i;
	
	i = _idlebw*(ie-_ebw)+1;
	xbwip = _bwspa(i)+ie*(_bwspb(i)+ie*(_bwspc(i)+ie*_bwspd(i)));
	
	return xbwip;
} //OK-2

/**********************************************************************
 * 
 *  Function: Linearly interpolated q(u;energy) surface      
 *	Input:	
 *      u -> angular variable                                      *
 *      ie -> -1/energy in eV^-1  --kinetic energy-- 
 *  Output:                                                      
 *      bw, broad screening parameter that gets flattest q surf
 *
 ***********************************************************************/
UREAL Cxq2Dip(UREAL u, UREAL ie, Txdpmq2D *xdpmq2D)
{
	UREAL ru,ou,rle,ole,ouole,xq2Dip;
	int32 i, j;
	
//	printf("u, iduq are %f, %f\n", u, iduq);
	ru = u*iduq;
	i = ru;
	ou = ru-i;
	i = i+1;
	rle = idleq*(ie-le0q);
	j = rle;
	ole = rle-j;
	j = j+1;
	ouole = ou*ole;
//	printf("ou, ole, ouole are %f, %f, %f\n", ou, ole, ouole);
	xq2Dip = q(i, j) * (1.0E0f - ou - ole + ouole) + q(i + 1, j) * (ou - ouole) + q(i, j + 1) * (ole - ouole) + q(i + 1, j + 1) * ouole;

	return xq2Dip;
} //OK-2

/**********************************************************************
 * 
 *  Function: Samples cos(theta) according to the G&S distribution. 
 *		Uses interpolated data for bw and the q surface and this
 *		latter quantity to perform a rejection procedure.      
 *	Input:	
 *      e -> kinetic energy in eV 
 *  Output:                                                      
 *      mu -> polar angle -cos(theta)-
 * 
 ***********************************************************************/
void Csamsca(UREAL e, UREAL *mu, Trxdpmbw *xdpmbw, Txdpmq2D *xdpmq2D,Trseed *rseed)
{
	UREAL bw, onebw, u, ie;
	UREAL tmp;

	ie = -1.0E0f / e;
	bw = Cxbwip(ie, xdpmbw);
	onebw = 1.0E0f + bw;
	//	printf("bw, onebw are %f, %f\n", bw, onebw);

	do
	{
		u = Crng(rseed);
		(*mu) = (onebw - u * (onebw + bw)) / (onebw - u);
		tmp = Cxq2Dip(u, ie, xdpmq2D);
		//printf("u, tmp are %f, %f\n", u, tmp);
	} while (Crng(rseed) > tmp);
	
}//OK-2

/**********************************************************************
 * 
 *  Function: Rotates a vector; the rotation is specified by giving
 *		the polar and azimuthal angles in the "self-frame", as
 *		determined by the vector to be rotated.       
 *	Input:	
 *      (u,v,w) -> input vector (=d) in the lab. frame
 *		costh -> cos(theta), angle between d before and after turn
 *      phi -> azimuthal angle (rad) turned by d in its self-frame 
 *  Output:                                                      
 *       (u,v,w) -> rotated vector components in the lab. frame
 * 
 ***********************************************************************/
void Crotate(UREAL *u, UREAL *v, UREAL *w, UREAL costh, UREAL phi)
{
	UREAL rho2,sinphi,cosphi,sthrho,urho,vrho,sinth,norm;
	
	rho2 = (*u)*(*u)+(*v)*(*v);
	norm = rho2+(*w)*(*w);
//     *** Check normalization:
	if(fabsf(norm-1.0E0f)>SZERO) //fabsf()
	{
		if(norm<ZERO)
		{
			printf("error:rotate:null vector cannot be renormalized");
			exit(1);
		}
//       *** Renormalize:
		norm = 1.0E0f/sqrtf(norm); //sqrtf()
		(*u) = (*u) * norm;
		(*v) = (*v) * norm;
		(*w) = (*w) * norm;	
	}
	
	sinphi = sinf(phi); //sinf()
	cosphi = cosf(phi); //cosf()
	
	//极端例子下可能出现1.0E0f-costh*costh<0;
	UREAL temp=1.0E0f-costh*costh;
	temp=(temp+abs(temp))/2.0E0f;
//     *** Case z' not= z:
	if(rho2 > ZERO) 
	{
		// sthrho = sqrtf((1.0E0f-costh*costh)/rho2); //sqrtf()
		sthrho = sqrtf((temp)/rho2); //sqrtf()
		urho =  (*u)*sthrho;
		vrho =  (*v)*sthrho;
		(*u) = (*u)*costh - vrho*sinphi + (*w)*urho*cosphi;
		(*v) = (*v)*costh + urho*sinphi + (*w)*vrho*cosphi;
		(*w) = (*w)*costh - rho2*sthrho*cosphi;	
	}else
	{
//     *** 2 especial cases when z'=z or z'=-z:
		// sinth = sqrtf(1.0E0f-costh*costh); //sqrtf()
		sinth = sqrtf(temp); //sqrtf()
		(*v) = sinth*sinphi;
		if ((*w)>0.0E0f)
		{
			(*u) = sinth*cosphi;
			(*w) = costh;
		} else
		{
			(*u) = -sinth*cosphi;
			(*w) = -costh;	
		}	
	}
} //OK-1

/**********************************************************************
 * 
 *  Function: Samples an energy loss according to Moller DCS      
 *	Input:	
 *      energy -> kinetic energy in eV  
 *  Output:                                                      
 *      eloss/energy, energy loss fraction
 * 
 ***********************************************************************/
/*
 * Hard to vectorization???
 */
UREAL Csammo(UREAL Lenergy, Tcutoff *cutoff,Trseed *rseed)
{
	UREAL ko,aelos2,kcut,a;
	UREAL sammo;
	
	kcut = _wcion / Lenergy;
	if (kcut < kcmax)
	{
		a = (1.0E0f - 1.0E0f / (1.0E0f + Lenergy * imc2));
		a = a * a;
		while (1 > 0)
		{
			if (Crng(rseed) * (1.0E0f + 2.5E0f * a * kcut) < 1.0E0f)
			{
				sammo = kcut / (1.0E0f - Crng(rseed) * (1.0E0f - 2.0E0f * kcut));
			}
			else
			{
				sammo = kcut + Crng(rseed) * (0.5E0f - kcut);
			}
			ko = sammo / (1.0E0f - sammo);
			aelos2 = a * sammo * sammo;
			if ((Crng(rseed) * (1.0E0f + 5.0E0f * aelos2)) < (1.0E0f + aelos2 + ko * (a + ko - 1.0E0f)))
			{
				break;
			}
		}
	}
	else
	{
		sammo = 0.0E0f;
	}
	return sammo;
} //OK-2

/**********************************************************************
 * 
 *  Function: gives the outgoing angle of the delta ray after a Moller
 *    interaction
 *	Input:	
 *      eloss -> energy loss --also kinetic energy of the recoil e-
 *		energy -> kinetic energy in eV
 *	Output:
 *		cos(theta') in the lab frame
 * 
 ***********************************************************************/
UREAL Csecmo(UREAL eloss, UREAL Lenergy)
{
	UREAL secmo;
	//sqrtf()
	secmo = sqrtf(eloss*(Lenergy+twomc2)/(Lenergy*(eloss+twomc2))); 
	
	return secmo;
} //OK-2


/**********************************************************************
 * 
 *  Function: Creates a new secondary electron from a Moller 
 *		interaction and stores its state in the secondary stack       
 *	Input:	
 *      elost -> kinetic energy of the secondary electron (eV)
 * 
 ***********************************************************************/
void Cputmol(UREAL elost, Tdpmpart *dpmpart, Trdpmstck *dpmstck,Trseed *rseed)
{
	if(_nsec>=maxns)
	{
		printf("putmol:error: Stack is full, enlarge cutoffs or maxns");
		exit(1);
	}
	_nsec = _nsec+1;
	_sptype(_nsec) = -1;
	_senerg(_nsec) = elost;
	_svx(_nsec) = vx;
	_svy(_nsec) = vy;
	_svz(_nsec) = vz;
	_sx(_nsec) = x;
	_sy(_nsec) = y;
	_sz(_nsec) = z;
	_sxvox(_nsec) = xvox;
	_syvox(_nsec) = yvox;
	_szvox(_nsec) = zvox;
	_sabsvx(_nsec) = absvox;
	Crotate(&_svx(_nsec), &_svy(_nsec), &_svz(_nsec), Csecmo(elost,energy), 
	Crng(rseed)*twopi);
} //OK-2


/**********************************************************************
 * 
 *  Function: Canibalized from PENELOPE too. 
 *		screening functions f1(b) and f2(b) in the bethe-heitler
 *		differential cross section for bremsstrahlung emission.
      
 *	Input:	
 *      elost -> kinetic energy of the secondary electron (eV)
 * 
 ***********************************************************************/
void Cschiff(UREAL *b, UREAL *f1, UREAL *f2)
{
	UREAL b2, a0;

	b2 = (*b) * (*b);
	*f1 = 2.0E0f - 2.0E0f * log(1.0E0f + b2); // logf()
	*f2 = (*f1) - 6.666666666666666E-1f;
	if (*b < 1.0E-10f)
	{
		(*f1) = (*f1) - twopi * (*b);
	}
	else
	{
		a0 = 4.0E0f * (*b) * atan2(1.0E0f, (*b)); // atan2f()
		(*f1) = (*f1) - a0;
		(*f2) = (*f2) + 2.0E0f * b2 * (4.0E0f - a0 - 3.0E0f * log((1.0E0f + b2) / b2)); // logf()
	}
	(*f2) = 0.5E0f * (3.0E0f * (*f1) - (*f2));
} //OK

/**********************************************************************
 * 
 *  Function: New routine to replace sambre in dpm  for more accurate
 *		sampling of Bremsstrahlung energies.
 *		Canibalized from PENELOPE.        
 *	Input:	
 *      elost -> kinetic energy of the secondary electron (eV)
 * 
 ***********************************************************************/
UREAL Csbrems(UREAL e, int32 m, Tcutoff *cutoff, Txbre *xbre,Trseed *rseed)
{
	//Random sampling of hard bremsstrahlung emission.
	int32 ic;
	UREAL toten, gam, em, ed, ec, ec1, ec2, emecs, xx, f2le;
	UREAL f00, f1, f2, f1c, f2c, f1d, f2d, ed1, bec, bed;
	UREAL eps, eps1, pa, pa0, pa1, pa2, pb, pb2;
	UREAL sbrems;
	
	UREAL xx1, xx2;
	
	toten=e+rev;
	gam=toten/rev;
	pa=bcb(m)*gam;
	em=e/toten;
	
	if(em > 9.9999E-1f) em=9.9999E-1f;

	//  ****  rejection functions.
	ed = (e - 5.0E0f * rev) / toten;
	ec = _wcbre / toten;
	ec2 = ec * ec;
	emecs = em * em - ec2;

	//  ****  low energy correction.
	xx = e / rev;
	xx1 = 1 + xx;
	xx1 = xx1 * xx1;
	xx2 = xx1 * (1 + xx);

	f2le = (4.650E0f - f0(m, 3) * (6.005E0f - f0(m, 3) * 2.946E0f)) / xx 
	- (3.242E1f - f0(m, 3) * (6.708E1f + f0(m, 3) * 3.906E0f)) / xx1 
	+ (2.033E1f + f0(m, 3) * (2.338E1f - f0(m, 3) * 7.742E1f)) / xx2;

	if (ec < ed)
	{
		ic = 2;
	}
	else
	{
		ic = 1;
	}

	f00 = f0(m, ic) + f2le;
	ec1 = 1.0E0f - ec;
	bec = ec / (pa * ec1);

	Cschiff(&bec, &f1, &f2);

	f1c = f1 + f00;
	f2c = (f2 + f00) * ec1;

	if (ec < ed)
	{
		f00 = f0(m, 1) + f2le;
		ed1 = 1.0E0f - ed;
		bed = ed / (pa * ed1);
		Cschiff(&bed, &f1, &f2);

		f1d = f1 + f00;
		f2d = (f2 + f00) * ec1;
		if (f1d > f1c)
			f1c = f1d;
		if (f2d > f2c)
			f2c = f2d;
	}

	pa1 = emecs * f1c;
	pa2 = 2.66666666666666E0f * logf(em / ec) * f2c; // logf()

LOOP1:
	if (Crng(rseed) * (pa1 + pa2) > pa1)
		goto LOOP2;
	eps = sqrtf(ec2 + Crng(rseed) * emecs); // sqrtf()
	pb = eps / (pa * (1.0E0f - eps));
	f1 = 2.0E0f - 2.0E0f * logf(1.0E0f + pb * pb); // logf()
	if (pb < 1.0e-8f)
	{
		f1 = f1 - pb * twopi;
	}
	else
	{
		f1 = f1 - 4.0E0f * pb * atan2f(1.0E0f, pb); // atan2f()
	}
	if (eps < ed)
	{
		f00 = f0(m, 2) + f2le;
	}
	else
	{
		f00 = f0(m, 1) + f2le;
	}
	if (Crng(rseed) * f1c > f1 + f00)
		goto LOOP1;
	goto LOOP3;

LOOP2:
	eps = powf((em / ec), Crng(rseed)); //
	eps = ec * eps;
	eps1 = 1.0E0f - eps;
	pb = eps / (pa * eps1);
	pb2 = pb * pb;
	f1 = 2.0E0f - 2.0E0f * logf(1.0E0f + pb2); // logf()
	f2 = f1 - 6.666666666666666E-1f;
	if (pb < 1.0E-8f)
	{
		f1 = f1 - pb * twopi;
	}
	else
	{
		pa0 = 4.0E0f * pb * atan2f(1.0E0f, pb); // atan2f()
		f1 = f1 - pa0;
		f2 = f2 + 2.0E0f * pb2 * (4.0E0f - pa0 - 3.0E0f * logf((1.0E0f + pb2) / pb2)); // logf()
	}
	f2 = 0.5E0f * (3.0E0f * f1 - f2);
	if (eps < ed)
	{
		f00 = f0(m, 2) + f2le;
	}
	else
	{
		f00 = f0(m, 1) + f2le;
	}

	if (Crng(rseed) * f2c > eps1 * (f2 + f00))
		goto LOOP1;

LOOP3:
	sbrems = eps * toten;

	return sbrems;
} // OK-2

/**********************************************************************
 * 
 *  Function: Creates a new secondary bremsstrahlung photon and 
 *		stores its state in the secondary stack        
 *	Input:	
 *      elost -> energy of the secondary photon (eV)
 * 
 ***********************************************************************/
void Cputbre(UREAL elost, Tdpmpart *dpmpart, Trdpmstck *dpmstck,Trseed *rseed)
{
	UREAL angsq2;
	_nsec = _nsec+1;
	if (_nsec>maxns)
	{
		printf("putbre:error: Stack is full, enlarge cutoffs or maxns");
	}
	_sptype(_nsec) = 0;
	_senerg(_nsec) = elost;
	_svx(_nsec) = vx;
	_svy(_nsec) = vy;
	_svz(_nsec) = vz;
	_sx(_nsec) = x;
	_sy(_nsec) = y;
	_sz(_nsec) = z;
	_sxvox(_nsec) = xvox;
	_syvox(_nsec) = yvox;
	_szvox(_nsec) = zvox;
	_sabsvx(_nsec) = absvox;
//		*** Polar angle set to mean value, i.e. no angular distribution used,
//		* besides, small angle approx is used for cosine:
	angsq2 = mc2sq2/(energy+mc2);
	Crotate(&_svx(_nsec),&_svy(_nsec),&_svz(_nsec),1.0E0f-angsq2*angsq2,
		Crng(rseed)*twopi);
} //OK-2

/**********************************************************************
 * 
 *  Function: Transports an electron until it either escapes from the
 *		universe or its energy drops below Eabs   
 *	Input:	
 *      electron initial state 
 *    Output:                                                      
 *      deposits energy in counters and updates secondary stack
 * 
 ***********************************************************************/
void Celectr(Tdpmpart *dpmpart, Tdpmvox *dpmvox, Tdpmmat *dpmmat, 
	Tdpmjmp *dpmjmp, Tdpmsrc *dpmsrc, Tdpmsub *dpmsub, dTdpmsim *dpmsim, 
	Trdpmscsr *dpmscsr, Trdpmlam *dpmlam, Tcgo *cgo,Tcgeom3 *cgeom3, 
	Trdpmstpw *dpmstpw, Trdpmscpw *dpmscpw, Trdpmlbr *dpmlbr, 
	Tdpmesc_new *dpmesc_new, Tdpmrpt *dpmrpt, Trxdpmbw *xdpmbw, Txdpmq2D *xdpmq2D,
	Trdpmstck *dpmstck, Tcutoff *cutoff, Txbre *xbre,Trseed *rseed)
{	
	UREAL inve, de, costhe;

	//*** Loading fuel before take off:
	ebefor = energy;
	fuelel = Cstepip(energy, dpmscsr);
	fuelxt = fuelel*Crng(rseed);
	fuelel = fuelel-fuelxt;
	inve = -1.0E0f/energy;
	fuelmo = -Clamoip(inve, dpmlam)*logf(Crng(rseed)); //logf
	fuelbr = -logf(Crng(rseed)); //logf

	//Set switches to steer the flight:
	modeel = 1;
	event = -1;

	//Repeat for every flight stop:
	while (1 > 0)
	{
		Cflight(dpmpart, dpmvox, dpmmat, dpmjmp, dpmsrc, dpmsub, cgo,
				cgeom3, dpmstpw, dpmscpw, dpmlbr, dpmesc_new, dpmrpt);
		//printf("It's fueling, event: %12d %12d\n", event, modeel);
		//printf("ISEED1,ISEED2,ICNT: %12d %12d %12d\n", ISEED1,ISEED2,ICNT);

		if (event == 2)
		{
#if DBG_OUT_ON == 1
			//Elastic stop for refueling:
			if (modeel == 0)
			{
				//End of second scattering substep, no real interaction:
				ebefor = energy;
				fuelel = Cstepip(energy, dpmscsr);
				fuelxt = fuelel * Crng(rseed);
				fuelel = fuelel - fuelxt;
				modeel = 1;
				event = 20;
			}
			else
			{
				//Elastic scattering event:
				//opt-qSingleMat ON: activate next line and deact next+1:
				//printf("Before samsca: %12d %12d %12d\n", ISEED1,ISEED2,ICNT);
				//printf("Before samsca: ebefor = %f \n", ebefor);
				Csamsca(ebefor, &costhe, xdpmbw, xdpmq2D, rseed);
				//printf("After samsca: %12d %12d %12d\n", ISEED1,ISEED2,ICNT);
				Crotate(&vx, &vy, &vz, costhe, twopi * Crng(rseed));
				fuelel = fuelxt;
				modeel = 0;
			}
#endif
		}
		else if (event == 3)
		{
#if DBG_OUT_ON == 1
			//Moller stop for refueling; note that:
			//no Moller is simulated below energy=2*cutoff;
			//energy/2 > de > cutoff(moller) >= Eabs  holds, then
			//EnergyFinal = energy-de can never be below Eabs;
			//sammo() can also return de=0.
			de = energy * Csammo(energy, cutoff, rseed);
			if (de > 0.0E0f)
			{
				Cputmol(de, dpmpart, dpmstck, rseed);
				energy = energy - de;
			}
			fuelmo = -Clamoip(-1.0E0f / energy, dpmlam) * logf(Crng(rseed)); // logf()
#endif
		}
		else if (event == 4)
		{
#if DBG_OUT_ON == 1
			//Bremsstrahlung stop for refueling; note that
			//de > cutoff(brem) >= Eabsph  OR  de=0 --see sambre():
			de = Csbrems(energy, _mat(absvox), cutoff, xbre,rseed);
			if (de > 0.0E0f)
			{
				Cputbre(de, dpmpart, dpmstck,rseed);
				energy = energy - de;
			}
			if (energy < _eabs)
			{
				//Determine if subEabs transport is needed:
				if (voxden > _subden)
				{
					Cscore(energy, dpmpart, dpmesc_new, dpmrpt);
				}
				else
				{
					Csubabs(dpmpart, dpmvox, dpmsub, dpmesc_new, cgo, dpmrpt, cgeom3);
				}
				break;
			}
			fuelbr = -logf(Crng(rseed)); // logf()
#endif
		}
		else
		{
			//*** Other event values indicate that history ended
			break;
		}
	}
} // OK

/**********************************************************************
 * 
 *  Function: Creates secondaries for annhilation photons, stores 
 *		states in the secondary stack    
 *	Input:	
 * 
 ***********************************************************************/
void Cputann(Tdpmpart *dpmpart, Trdpmstck *dpmstck, Tdpmsrc *dpmsrc,Trseed *rseed)
{
	UREAL sinthe, phi;
	int32 j;
	
	vz = 2.0E0f*Crng(rseed)-1.0E0f;
	sinthe = sqrtf(1.0E0f-vz*vz); //sqrtf()
	phi = twopi*Crng(rseed); 
	vy = sinthe*sinf(phi); //sinf()
	vx = sinthe*cosf(phi); //cosf()
	for(j=1; j<=2; j++)
	{
		_nsec = _nsec+1;
		if (_nsec>maxns)
		{
			printf("putann:error: stack is full, enlarge maxns.");
			exit(1);
		}
		_sptype(_nsec) = 0;
		_senerg(_nsec) = mc2;
		_svx(_nsec) = vx;
		_svy(_nsec) = vy;
		_svz(_nsec) = vz;
		_sx(_nsec) = x;
		_sy(_nsec) = y;
		_sz(_nsec) = z;
		_sxvox(_nsec) = xvox;
		_syvox(_nsec) = yvox;
		_szvox(_nsec) = zvox;
		_sabsvx(_nsec) = absvox;
	}
	_svx(_nsec) = -_svx(_nsec);
	_svy(_nsec) = -_svy(_nsec);
	_svz(_nsec) = -_svz(_nsec);
} //OK-1

/**********************************************************************
 * 
 *  Function: Mean free path prepared to play the Woodcock trick  
 *	Input:
 *		e -> kinetic energy in eV
 *	Output:
 *		Minimum mean free path in cm
 *
 ***********************************************************************/
UREAL Clamwck(UREAL e, Trdpmwck *dpmwck)
{
	int32 i;
	UREAL lamwck;
	
	i = (e-_wcke0)*_idlewk+1;
	lamwck = _a0wck(i)+_a1wck(i)*e;
	
	return lamwck;
}//Ok-1

/**********************************************************************
 * 
 *  Function: Inverse Compton mean free path --3spline interpolation  
 *	Input:
 *		matid -> material id# 
 *		e -> kinetic energy in eV
 *	Output:
 *		Inverse total mean free path in cm^2/g 
 *
 ***********************************************************************/
UREAL Cicptip(int32 Lmatid, UREAL e, Trdpmcmp *dpmcmp)
{
	int32 i;
	UREAL icptip;
	
	i = _idlecp*(e-_ecmpt)+1;
	icptip = _compta(Lmatid,i) + e*(_comptb(Lmatid,i) + e*(_comptc(Lmatid,i) 
		+ e*_comptd(Lmatid,i)));
		
	return icptip;
} //OK-1

/**********************************************************************
 * 
 *  Function: Photon total inverse mean free path --3spline 
 *		interpolation  
 *	Input:
 *		matid -> material id#
 *		e -> kinetic energy in eV
 *	Output:
 *		Total inverse mean free path in cm^2/g
 *
 ***********************************************************************/
UREAL Citphip(int32 Lmatid, UREAL e, Trdpmlph *dpmlph)
{
	int32 i;
	UREAL itphip;
	
	i = _idleph*(e-_elaph)+1;
	itphip = _lampha(Lmatid,i) + e*(_lamphb(Lmatid,i) + e*(_lamphc(Lmatid,i)
		+e*_lamphd(Lmatid,i)));
	return itphip;
} //Ok-1

/**********************************************************************
 * 
 *  Function: Compton angular deviation of the secondary electron.  
 *	Input:
 *		energy -> photon energy in eV
 *		efrac -> fraction of initial energy kept by 2nd photon
 *		costhe-> photon scattering angle
 *  Output:                                                      *
 *  	-> cos(theta) of the 2nd electron  
 *
 ***********************************************************************/
UREAL Ccomele(UREAL Lenergy, UREAL efrac, UREAL costhe)
{
	UREAL e0, comele;
	
	e0 = Lenergy*imc2;
	/* sqrtf() */
	comele = (1.0E0f + e0) * sqrt((1.0E0f - efrac) / (e0 * (2.0E0f + e0 * (1.0E0f - efrac))));
		
	return comele;
} //OK-1

/**********************************************************************
 * 
 *  Function: Creates a secondary electron from a Compton interaction 
 *		and stores its state in the secondary stack.  
 *	Input:
 *		elost -> energy of the secondary electron being created
 *		efrac -> fraction of initial energy kept by 2nd photon
 *		phi -> athimutal angle
 *		costhe -> photon scattering angle
 *
 ***********************************************************************/
void Cputcom(UREAL elost, UREAL efrac, UREAL phi, UREAL costhe, 
	Tdpmpart *dpmpart, Trdpmstck *dpmstck)
{
	_nsec = _nsec+1;	
	if (_nsec>maxns)
	{
		printf("putcom:error: Stack is full, enlarge cutoffs or maxns.");
		exit(1);
	}
	_sptype(_nsec) = -1;
	_senerg(_nsec) = elost;
	_svx(_nsec) = vx;
	_svy(_nsec) = vy;
	_svz(_nsec) = vz;
	_sx(_nsec) = x;
	_sy(_nsec) = y;
	_sz(_nsec) = z;
	_sxvox(_nsec) = xvox;
	_syvox(_nsec) = yvox;
	_szvox(_nsec) = zvox;
	_sabsvx(_nsec) = absvox;
	
	Crotate(&_svx(_nsec),&_svy(_nsec),&_svz(_nsec),
		Ccomele(energy,efrac,costhe), phi);
} //OK-1

/**********************************************************************
 * 
 *  Function: Inverse pair production mean free path --3spline interpol  
 *	Input:
 *		matid -> material id#
 *		e -> kinetic energy in eV
 *  Output:
 *  	Inverse total mean free path in cm^2/g 
 *
 ***********************************************************************/
/*
 * Change pairpa,b,c,d的存储方式,连续存储？？？
 */
UREAL Cipapip(int32 Lmatid, UREAL e, Trdpmpap *dpmpap)
{
	int32 i;
	UREAL ipapip;
	
	if (e>_epair)
	{
        i = _idlepp*(e-_epair)+1;
        ipapip = _pairpa(Lmatid,i)+e*(_pairpb(Lmatid,i) + 
			e*(_pairpc(Lmatid,i)+e*_pairpd(Lmatid,i)));
	} else
	{
		ipapip = 0.0E0f;
	}
		
	return ipapip;
} //OK-1

/**********************************************************************
 * 
 *  Function: Samples a Compton event following Klein-Nishina DCS  
 *	Input:
 *		energy -> photon energy in eV
 *  Output:
 *  	efrac -> fraction of initial energy kept by 2nd photon
 *		costhe -> cos(theta) of the 2nd photon 
 *
 ***********************************************************************/
void Ccomsam(UREAL Lenergy, UREAL *efrac, UREAL *costhe,Trseed *rseed)
{
	UREAL e0,twoe,twoe1,kmin2,loge,mess;
	UREAL a1, b1;
	
	e0 = Lenergy*imc2;
	twoe = 2.0E0f*e0;
/* 	twoe1 = 1.0E0f+twoe;
	kmin2 = 1.0E0f/(twoe1*twoe1);
	loge = logf(twoe1); //logf() */
	kmin2 = 1.0E0f/powf((1.0E0f+twoe),2);  //powf
	loge = logf(1.0E0f+twoe);            
	
	do {
		a1 = Crng(rseed)*(loge+twoe*(1.0E0f+e0)*kmin2);
		if(a1 < loge)
		{
			*efrac= expf(-Crng(rseed)*loge);//expf
		} else
		{
			*efrac = sqrtf(kmin2+Crng(rseed)*(1.0E0f-kmin2));
		}
		mess = e0*e0*(*efrac)*(1.0E0f+(*efrac)*(*efrac));
		a1 = Crng(rseed)*mess;
		b1 = mess-(1.0E0f-(*efrac))*((1.0E0f+twoe)*(*efrac)-1.0E0f);
	}while(a1 > b1);
		
	*costhe = 1.0E0f-(1.0E0f-(*efrac))/((*efrac)*e0);
} //OK-1

/**********************************************************************
 * 
 *  Function: Transports a photon until it either escapes from the  
 *		universe or its energy drops below EabsPhoton  
 *	Input:
 *		photon initial state
 *	Output:
 *		deposits energy in counters and updates secondary stack
 * 
 ***********************************************************************/
void Cphoton(Tdpmpart *dpmpart, Trdpmstck *dpmstck, Tdpmvox *dpmvox, 
	Tdpmsrc *dpmsrc, Tdpmesc_new *dpmesc_new, Trdpmwck *dpmwck, Trdpmcmp *dpmcmp,
	Tdpmrpt *dpmrpt, Tcgeom3 *cgeom3, Trdpmlph *dpmlph, Trdpmpap *dpmpap,Trseed *rseed)
{
	UREAL s, lammin, lamden, prob, randno, costhe, efrac, de;
	UREAL phi;
	UREAL epair[2];
	int32 j;
	
//     *** Loop until it either escapes or is absorbed:
	while(1>0){
//     *** Get lambda from the minimum lambda at the current energy:
		lammin = Clamwck(energy, dpmwck);
		s = -lammin*logf(Crng(rseed)); //logf()
		x = x+s*vx;
		y = y+s*vy;
		z = z+s*vz;
		Cwhere(dpmpart, cgeom3, dpmvox);
#if DBG_OUT_ON == 0			
		printf("zvox, yvox, xvox are %12d %12d %12d \n", zvox, yvox, xvox);
#endif
		if(absvox==0)
		{
			break;
		}
//       *** Apply Woodcock trick:
        lamden = lammin*_dens(absvox);
        prob = 1.0E0f-lamden*Citphip(_mat(absvox),energy,dpmlph);
        randno = Crng(rseed);
		
//       *** No real event; continue jumping:
        if (randno < prob) continue;

//       *** Compton:
        prob = prob+lamden*Cicptip(_mat(absvox),energy,dpmcmp);
		if (randno<prob)
		{
//c         *** opt-IncohScat-ON ->  Activate sincoh() and deact comsam():
//c         * call sincoh(_mat(absvox),energy,efrac,costhe)
#if DBG_OUT_ON == 0	
			cnt1=cnt1+1;
#endif
			Ccomsam(energy,&efrac,&costhe,rseed);
			de = energy*(1.0E0f-efrac);
			phi = twopi*Crng(rseed);
			if (de<_eabs)
			{
				Cscore(de, dpmpart, dpmesc_new, dpmrpt);
#if DBG_OUT_ON == 0	
			cnt1=cnt1+1;
#endif				
			} else
			{
//           *** Create a secondary electron in the stack:				
				Cputcom(de,efrac,phi+pi,costhe,dpmpart,dpmstck);		
			}
			energy = energy-de;
			if(energy < _eabsph)
			{
				Cscore(energy, dpmpart, dpmesc_new, dpmrpt);
#if DBG_OUT_ON == 0	
			cnt2=cnt2+1;
#endif				
				break;
			}
			Crotate(&vx,&vy,&vz,costhe,phi);
			continue;
		}

//       *** Pair production:		
		prob = prob+lamden*Cipapip(_mat(absvox),energy,dpmpap);
		if(randno < prob)
		{
#if DBG_OUT_ON == 0	
			cnt2=cnt2+1;
#endif			
			epair[0] = Crng(rseed)*(energy-twomc2);
			epair[1] = energy-twomc2-epair[0];
			for(j=1; j<=2; j++)
			{
				if(epair[j-1]>_eabs)
				{
					_nsec=_nsec+1;
					if (_nsec>maxns)
					{
						printf("photon:error: Stack is full, enlarge maxns.");
						exit(1);
					}
					_senerg(_nsec) = epair[j-1];
					_svx(_nsec) = vx;
					_svy(_nsec) = vy;
					_svz(_nsec) = vz;
					_sx(_nsec) = x;
					_sy(_nsec) = y;
					_sz(_nsec) = z;
					_sxvox(_nsec) = xvox;
					_syvox(_nsec) = yvox;
					_szvox(_nsec) = zvox;
					_sabsvx(_nsec) = absvox;
					if(j==1)
					{
						_sptype(_nsec) = -1;
					} else
					{
						_sptype(_nsec) = +1;
					}
				} else
				{
					Cscore(epair[j-1], dpmpart, dpmesc_new, dpmrpt);
#if DBG_OUT_ON == 0	
			cnt3=cnt3+1;
#endif						
					if(j==2)
					{
						Cputann(dpmpart, dpmstck, dpmsrc,rseed);
					}
				}
			}
			break;
		}
//    *** In any other case, photoelectric absorption occurs:
		Cscore(energy, dpmpart, dpmesc_new, dpmrpt);
#if DBG_OUT_ON == 0	
			cnt3=cnt3+1;
#endif			
		break;
	}
	
} //OK-1

/**********************************************************************
 * 
 *  Function: Retrieves a particle from the secondary stack and fills 
 *		up the current particle common with its state   
 *  Output:
 *  	-> 0 if there was no particle in the stack, 1 else  
 *
 ***********************************************************************/
int32 Cscndry(Tdpmpart *dpmpart, Trdpmstck *dpmstck)
{
	int32 scndry;
	
	if(_nsec==0)
	{
		scndry = 0;
	}else
	{
		scndry = 1;
		energy = _senerg(_nsec);
		vx = _svx(_nsec);
		vy = _svy(_nsec);
		vz = _svz(_nsec);
		x = _sx(_nsec);
		y = _sy(_nsec);
		z = _sz(_nsec);
		xvox = _sxvox(_nsec);
		yvox = _syvox(_nsec);
		zvox = _szvox(_nsec);
		absvox = _sabsvx(_nsec);
		ptype = _sptype(_nsec);
		_nsec = _nsec-1;
	}

	return scndry;
} //OK-1
/**********************************************************************
 * 
 *  Function: Dumps all tmp counters into the corresponding final 
 *		counters   
 *  Output: 
 *
 ***********************************************************************/
void Cdumpe(Tdpmesc *dpmesc)
{
	int32 i;
	
	for(i=1; i<=nxyz; i++)
	{
		escore(i) = escore(i)+etmp(i);
		escor2(i) = escor2(i)+etmp(i)*etmp(i);
		etmp(i) = 0.0E0f;
	}
}

void Cdumpe_new(Trdpmesc *dpmesc,Tdpmesc_new *dpmesc_new)
{
	int32 i;

	for(i=0; i<tail; i++)
	{
		_escore(esidx(i)) = _escore(esidx(i))+estmp(i);
		_escor2(esidx(i)) = _escor2(esidx(i))+estmp(i)*estmp(i);
	}
}
/**********************************************************************
 * 
 *  Function: 计算每个进程分配的粒子数 
 *		   
 *  Output: 
 *
 ***********************************************************************/
int32 Ccomp_np(int32 nprtcl_all,int32 gloabl_id,int32 gloabl_num)
{
  int32 nprtcl;
  long long tp_np=nprtcl_all;
  long long tp_gid=gloabl_id;
  long long tp_gnm=gloabl_num;
  
  // 根据ID号计算分配得到的粒子数
  nprtcl = int32 (tp_np * (tp_gid+1) / tp_gnm  - tp_np * (tp_gid) / tp_gnm);

  return(nprtcl);
}
