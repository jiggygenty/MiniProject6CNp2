#include<stdio.h>
#include<stdlib.h>
#include<math.h>

double u0(double n,double x)
{
    return(cos((n-0.5)*x));
}
double uexact(double n,double t,double x)
{
    return(exp(-pow(n-0.5,2)*t)*cos((n-0.5)*x));
}

void main()
{
double n=1;
double D=1.0;
double F=0.0;
double dx=M_PI/20;
double xi=0.0;
double xf=M_PI;
double tMAX=2.0;
double lambda=0.4;
double dt=lambda*dx*dx/D;
int i,j;
int N=(int)ceil((xf-xi)/dx);

int M=(int)ceil(tMAX/dt);
double u[M+1][N+1];
double a1,a2;
double b1,b2;
double c1,c2;
double arrayres[N+1];
double alpha[N+1];
double g[N+1];
double exact[M+1][N+1];
double error[M+1][N+1];

FILE *fptr;
FILE *fr;
fptr=fopen("CNdata4p2.csv","w");
fr=fopen("CNerror4p2.csv","w");

for(i=0;i<=M;i++)
{
    for(j=0;j<=N;j++)
    {
    printf("%lf ",uexact(n,(i*dt),(j*dx)));
    exact[i][j]=uexact(n,(i*dt),(j*dx));
    }
    printf("\n");
}
printf("\n");
printf("\n");

for(j=0;j<=N;j++)
{
    u[0][j]=u0(n,j*dx);
}

a1=(1-lambda);
c1=lambda/2.0;
b1=lambda/2.0;

a2=1+lambda;
c2=-lambda/2.0;
b2=-lambda/2.0;

for(i=0;i<=M;i++)
{
    u[i][N]=0;
}

for(i=1;i<=M;i++)
{
alpha[0]=a2;
g[0]=(c1)*u[i-1][1]+((a1)*u[i-1][0]+(b1)*u[i-1][1]);


for(j=1;j<=N-1;j++)
{
        if(j==1)
        {
        alpha[j]=a2-2*(b2*c2)/alpha[j-1];
        }
        else
        alpha[j]=a2-(b2*c2)/alpha[j-1];
        g[j]=(c1)*u[i-1][j-1]+((a1)*u[i-1][j]+(b1)*u[i-1][j+1])-b2/alpha[j-1]*g[j-1];
}

u[i][N-1]=g[N-1]/alpha[N-1];
for(j=N-2;j>=1;j--)
{
    u[i][j]=(g[j]-(c2*u[i][j+1]))/alpha[j];
}
u[i][0]=(g[0]-2*c2*u[i][1])/alpha[0];

}

for(i=0;i<=M;i++)
{
    for(j=1;j<=N-1;j++)
    {
      error[i][j]=100.0*(fabs(u[i][j]-exact[i][j]))/exact[i][j];
      fprintf(fr,"%lf,",error[i][j]);
    }
    fprintf(fr,"\n");
}
for(i=0;i<=M;i++)
{
for(j=0;j<=N;j++)
{
    printf("%lf ",u[i][j]);
}
printf("\n");
}
for(i=0;i<=M;i++)
{
    for(j=0;j<=N;j++)
    {
        fprintf(fptr,"%lf,",u[i][j]);
    }
    fprintf(fptr,"\n");
}
fclose(fptr);
fclose(fr);
}
