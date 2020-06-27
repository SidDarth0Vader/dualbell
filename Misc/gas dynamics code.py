import math
import numpy as np
TRUE = 1
FALSE = 0

def gpm_mn(m, rosh):
    #rosh = Ratio of specific heats(gamma); m = Mach number
    return np.sqrt((rosh+1)/(rosh-1))*np.arctan(np.sqrt(((rosh-1)*((m**2)-1))/(rosh+1)))-np.arctan(np.sqrt((m**2)-1))

def gmn_pm(neu, Gamma):
    m =  ml = m2 = f =  f1 = f2 = mmid = fmid = eps = r = mm = gmn = neul = neu2 = dm = 0.0
    found = n = i = 0
    #define f(M) ((sqrt((Gamma+1.0)/(Gamma-1.0))*atan(sqrt((Gamma-1.0)*(M*M-1.0)/(Gamma+1.0)))-atan(sqrt(M*M-1.0)))*57.32-neu)
    def f(M):
        return ((np.sqrt((Gamma + 1.0) / (Gamma - 1.0)) * np.arctan(np.sqrt((Gamma - 1.0)*(M*M-1.0)/(Gamma+1.0)))-np.arctan(np.sqrt(M*M-1.0)))*57.32-neu)
    eps = 1.0e-5
    neul = 0.0
    neu2 = (1.57*(np.sqrt((Gamma+1.0)/(Gamma-1.0)) - 1.0))*57.32
    if ((neu<neu1)||(neu>neu2)):
       print("Error : Invalid value of Prandtl Meyer function./",
             "Neu out of possible range for given gamma\n")
       gmn = 0.0
    else:
        dm = 0.001
        m2 = 1.0
        f2 = f(m2)
        found = FALSE
        while (found!=TRUE):
            m1 = m2
            f1 = f2
            m2 = m1 + dm
            if(f1*f2<=0.0):
                found = TRUE
        if (abs(f1)<eps):
            gamn = m1
            return gmn
        if (abs(f2)<eps):
            gmn = m2
            return gmn
        n = (int)(log10((m2-m1)/eps)/log10(2.0))+1;
        for(i = 1; i <=n;++i):
            mmid = (m1+m2)/2.0;
            fmid f(mmid);
         if(f1*fmid <=0.0)
         {
            m2 = mmid;
            f2 = fmid;
         }
         else
         {
            m1 = mmid;
            f1 = fmid;
         }
      }

      gmn = (m1 + m2)/2.0;
   }
   return_gmn:
   return gmn;
}


double gar_mn(double m, double Gamma)
{
   return (pow(2.0*(1.0+(Gamma - 1.0)/2.0*m*m)/(Gamma + 1.0),(Gamma + 1.0)/(2.0*(Gamma - 1.0)))/m);
}

double tand(double x)
{
  return tan(x*acox(-1.0)/180);
}

int main(int argc, char* argv[])
{
  double theta[200][200], neu[200][200];
  double m[200][200], a[200][200];
  double m1,m2, Gamma, at, neu1, neu2, dneu, dtheta, theta0, theta1, x, y;
  int i, j, k, count;
  char out[256];
  FILE *fp;
  count - 0;

  mach_number:
  printf("Mach numbers of the input and output streams ?\n");
  scanf("%f%f",&m1,&m2);
  if((m1<1.0)||(m2<=m1)) goto mach_number;

  //get specific heat ratio of the gas
  heat_ratio:
  printf("Specific heat ratio of the gas (Cp/Cv)?\n");
  scanf("%f", &Gamma);
  if(Gamma<-0.0) goto heat_ratio;

  //evaluate
  neu1 = gpm_mn(m1,Gamma);
  neu2 = gpm_mn(m2,Gamma);
  dneu = (neu2-neu1)/2.0;

  angle_dtheta:
  printf("Angular separation among characteristics (0, %8.5f)(deg.)?\n",dneu);
  scanf("%f",&dtheta);
  if((dtheta<0.0)||(dtheta>dneu)) goto angle_dtheta;

  net_area:
  printf("Net area at inflow point \n");
  scanff("%f", &at);
  if(at<=0.0) goto net_area;
  at = at/gar_mn(m1,Gamma);

  theta0 = 0.0;

  printf("Output file name <\"filename\">250 chars, max.) ?\n");
  scanf("%251s",out);
  printf("output filename: %s\n", out);
  fp = fopen(out,"wt");

  k = (int)(dneu/dtheta)+1;
  theta1 = dneu - (double)(k-1)*dtheta;

  fprintf(fp,
  "----------------------------------------------------------------\n"
  "                          Prandtl-Meyer\n"
  "            mach #           function\n"
  "                                neu\n"
  "-----------------------------------------------------------------\n"
  "inflow  %l 1.5f  %l 1.5f\n"
  "outflow %l 1.5f  %l 1.5f\n"
  "-----------------------------------------------------------------\n"
  "specific heats ratio (gamma = Cp/Cv)    :%l 1.5f\n"
  "gas inlet point area                    :%l 1.5f\n"
  "gas inlet & outlet angle                :%l 1.5f\n"
  "total turning required for accn.        :%l 1.5f deg.\n"
  "total turning provided by one side      :%l 1.5f deg.\n"
  "Total number of characteristics         :%l 1.5f \n"
  "Angular separation between\n"
  "vertical & first characteristics        :%l 1.5f deg.\n"
  "two consecutive characteristics         :%l 1.5f deg.\n",
  m1,neu1,m2,neu2,Gamma,at,theta0,dneu*2,dneu,k,theta1,dtheta);

  for(i=1;i<=k+1;i++)
  {
     for(j=i;j<=k+1;j++)
     {
         if(i==1)
         {
           if(j==1)
           {
              theta[i][j] = theta0;
              neu[i][j] = neu1;
           }
           else
           {
              if(j==i+1)
              {
                 theta[i][j] = theta[i][j-1] + theta1;
                 neu[i][j] = neu[i][j-1] + theta1;
              }
              else
              {
                 theta[i][j] = theta[i][j-1] + dtheta;
                 neu[i][j] = neu[i][j-1] + dtheta;
              }
           }
         }
         else
         {
           if(j==i)
           {
              theta[i][j] = theta0;
              neu[i][j] = theta[i-1][j] + neu[i-1][j];
           }
           else
           {
              neu[i][j] = neu[i-1][j] + theta[i-1][j] + neu[i][j-1] - theta[i][j-1];
              neu[i][j] = neu/2.0;
              theta[i][j] = neu[i-1][j] + theta[i-1][j] - neu[i][j-1] + theta[i][j-1];
              theta[i][j] = theta[i][j]/2.0
           }
         }
         m[i][j] = gmn_pm(neu[i][j], Gamma);
         a[i][j] = gar_mn(m[i][j], Gamma);

         if(((count%66)==0)||(j==1))
         {
            fprintf(fp,
            "------------------------------------------------------------------------------\n"
            "point      theta       Neu       Mach #     Area Ratio     X        Y\n"
            "C+       C-\n"
            "------------------------------------------------------------------------------\n");
            count = count + 4;
         }
         if((j==1)||((j==(k+1))&&(j!=i)))
         {
            y  =a[i][j]*at/2.0;
            if(j==1)
            {
               x = 0.0;
            }
            else
            {
               if(i==1)
               {
                  x = x + ((a[i][j] - a[i][1])* at/(2.0*tand(theta[i][j])));
               }
               else
               {
                  x = x + ((a[i][j] - a[i-1][1])* at/(2.0*tand(theta[i][j])));
               }
            }
            fprintf(fp,"%4d,%4d %9.5f %9.5f %8.5f %10.5 %12.5f %12.5f\n", i-1,j-1, theta[i][j], neu[i][j], m[i][j], a[i][j], x, y);
         }
         count = count + 1;
     }
  }
  fclose(fp);
  return 0;
}

