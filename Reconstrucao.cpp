#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define sqr(a) ((a)*(a))

//===================================================
// filtro espacial ram-lak
//===================================================
double h(int n, double dt) {
    // n: indice da posicao do filtro
    double resp; // valor do filtro em funcao da posicao

    if (n && n%2) resp = -1.0/sqr(n*M_PI*dt); // se impar
    else if (n) resp = 0;                     // se par
    else resp = 1/sqr(2*dt);                  // se nulo
    return resp;
}

//===================================================
// funcao convolucao entre projecao e ram-lak
//===================================================
double *Convolucao(int nn, double *p, double dt) {
    // nn: numero de raios-soma na projecao
    // *p: projecao
    // dt: resolucao espacial
    double *q=new double [nn]; // resultado da convolucao
    int     n, k;              // indices

    for ( n=0; n<nn; n++ ) {
        q[n] = 0;
        for ( k=0; k<nn; k++ ) {
            q[n] += p[k]*h(n-k,dt); // termo da convolucao
        }
    }
    return q;
}

int main(void) {
    FILE *fp;
    char fnome[64],str[3], reconst[64], inverte[64];
    int dr,dq;
    int  i, j, k;   // indices
    int  n,         // n raios-soma
         m,         // m projecoes
         t,         // t variavel auxiliar
         r,an;       //projecoes e raios-soma do arquivo fonte
    double **proj,  // projecoes de entrada
             dteta; // passo angular
    double **rec,   // matriz de reconstrucao
            *q,     // projecao filtrada
             val;   // maior valor reconstruido
    int ii, jj;
    
    dq=dr=1;
    r=600;
    an=600;

/*

//===================================================
//funcao de inversao de contagem para projecao
//===================================================
	printf("abrindo arquivo para inversao do PGM\n");
	sprintf(inverte,"inverte.pgm");
	fp = fopen(inverte,"rt");
	if (fp==NULL){
        printf ("erro na abertura do arquivo de inversao\n");
        exit(0);
    }
    fscanf(fp,"%s %d %d %lf\n",str,&n,&m,&val);
     
	proj = new double*[m];       // m projecoes
    for (i=0; i<m; i++ ) {
       proj[i] = new double[n]; // n raios-soma por projecao
    }
	
	printf("lendo as projecoes a serem invertidas\n");
    // leitura das projecoes
    for (i=0; i<m; i++ ) {
        for (j=0; j<n; j++ ) {
            fscanf(fp,"%lf",&proj[i][j]);
        }
    }
    fclose(fp);
    
    printf("salvando material invertido\n");
    // salvando material invertido
    fp = fopen("invertido.pgm","wt");
    	printf("arquivo invertido.pgm\n");
        fprintf(fp,"P2 %d %d %lf\n", n, m, val);
        for (i=0; i<m; i++ ) {
            for (j=0; j<n; j++){
            	
            	//fprintf(fp,"%4d",(int)((val)-(proj[i][j])));
            	
            	if(i<(0.05*m) || i>(m-(0.05*m))){
					fprintf(fp,"%4d",0);
				}else if(j<(0.05*n) || j>(n-(0.05*n))){
					fprintf(fp,"%4d",0);
				}else{
                	fprintf(fp,"%4d",(int)((val)-(proj[i][j])));
                }
            }
            fprintf(fp,"\n");
        }
    fclose(fp);
    
*/
	printf("abrindo sinograma para reconstrucao\n");
    sprintf(fnome,"sino.pgm");
    fp = fopen(fnome,"rt");
    if (fp==NULL){
        printf ("erro na abertura do sinograma\n");
        exit(0);
    }
       
    
        // str: "P2"; n: raios-soma; m: projecoes; val: maior valor (nao usado)
        fscanf(fp,"%s %d %d %lf\n",str,&n,&m,&val);

        // passo angular
        dteta = 2*M_PI/(m);

        // alocando espaco para as projecoes
        proj = new double*[m];       // m projecoes
        for (i=0; i<m; i++ ) {
            proj[i] = new double[n]; // n raios-soma por projecao
        }

        // leitura das projecoes
        for (i=0; i<m; i++ ) {
            for (j=0; j<n; j++ ) {
                fscanf(fp,"%lf",&proj[i][j]);
            }
        }
    fclose(fp);

    //=========================================
//    printf("raios-soma: %d\tprojecoes: %d\n",n,m);
//    for (i=0; i<m; i++ ) {
//        for (j=0; j<n; j++ ) {
//            printf("%lf ",proj[i][j]);
//        }
//        printf("\n");
//    }
    //=========================================
	
	printf("realizando reconstrucao\n");
    // alocando espaco para a reconstrucao
    rec = new double*[n];
    for (i=0; i<n; i++ ) {
        rec[i] = new double[n];
        for (j=0; j<n; j++ )
            rec[i][j] = 0;   // zerando a matriz
    }

    for ( i=0; i<m; i++ ) {
        // filtragem da projecao usando ram-lak
        q = Convolucao(n,proj[i],1);

        // retroprojecao da projecao filtrada (a alma do negocio)
        for (ii=0; ii<n; ii++ ) {
            for (jj=0; jj<n; jj++ ) {
                // para cada elemento da matriz de reconstrucao,
                // projeta a coordenada sobre a projecao
                t=(int)(n/2+(jj-n/2)*cos(i*dteta)+(ii-n/2)*sin(i*dteta));
                // verifica se a posicao eh valida e ...
                if (t>=0 && t<n)  // ... acumula (pincel)
                    rec[ii][jj] += q[t]*dteta;
            }
        }
        // descarta a projecao filtrada usada
        delete [] q;
    }

    // aqui acabou a reconstrucao; falta salvar o resultado...

    // busca o maior valor reconstruido
    val = 0;
    for (ii=0; ii<n; ii++ ) {
        for (jj=0; jj<n; jj++ ) {
            if (val<rec[ii][jj]) val = rec[ii][jj];
        }
    }
    
    
    // zera os valores inferiores a 10% de val (so pra limpar a imagem)
    for (ii=0; ii<n; ii++ ) {
        for (jj=0; jj<n; jj++ ) {
            if (rec[ii][jj]<0.1*val) rec[ii][jj] = 0;
        }
    }
    
    
    printf("salvando arquivo reconstruido\n");
    // salva a reconstrucao reescalando os valores para 255
    sprintf(reconst,"reconstrucao_pr%d_rs%d.pgm",(n),(m));
    fp = fopen(reconst,"wt");
        fprintf(fp,"P2 %d %d %d # %lf\n", n, m, 255, val);
        for (ii=0; ii<n; ii++ ) {
            for (jj=0; jj<n; jj++ ) {
                fprintf(fp,"%4d",(int)(255*rec[ii][jj]/val));
            }
            fprintf(fp,"\n");
        }
    fclose(fp);

    // acabou...
    for (i=0; i<m; i++ ) {
        delete [] proj[i];  // libera memoria
    }
    delete [] proj;   // libera memoria

    for (i=0; i<n; i++ ) {
        delete [] rec[i];  // e mais memoria
    }
    delete [] rec;  // e mais memoria

    return 0; // fim
}


