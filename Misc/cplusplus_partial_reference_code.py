import math
import sys

def main():
    float m, nu, theta[500], Nu[500], MinusK[500], PlusK[500], z[500], thetamax, t, a, rem=0.0, remain[2];
    int n, dt = 0, i, p, j=0, jo ,l=0 , k=2, pts[200], total=0, points[200], pnts[200], wallpts[200], intpts[200], y[200], first = 0, second = 1, x=2, q=1;
    print("Enter the Design Exit Mach Number.\n")
    scanf("%f",&m);
    printf("Enter the Flow Deflection Angle.\n");
    scanf("%f",&t);
    printf("Enter the Number of Characteristics to be considered.\n");
    scanf("%d",&n);
    rem = fmod(t,(n-1));
    printf("\nRem=%f\n",rem);
    a=t-rem;
    remain[0] = rem;
    printf("\nThe value of a=%f",a);
    dt = a/(n-1);
    printf("\nDelta Theta = %d\n\n", dt);
    theta[pts[0]] = rem;
    printf("\nThe points on the First Left Running Characteristic line are as follows:\n\n");
    for(i=0;i<n;i++)
    {  pts[i]= ++l;
       //printf("%d\n",pts[i]);
       printf("Pts[%d] = %d\n", i , pts[i]);
    }
    printf("\n Flow properties of the points on the First Left Running Characteristic line:\n\n");
    for(i=0; i<n; i++)
    {  theta[pts[i]] = rem + (dt*j);
       Nu[pts[i]] = theta[pts[i]];
       MinusK[pts[i]] = theta[pts[i]]+ Nu[pts[i]];
       PlusK[pts[i]] = theta[pts[i]] - Nu[pts[i]];
       j=j+1;
       printf("Point %d\tMinusK[%d] = %f\tPlusK[%d] = %f\tTheta[%d] = %f\tNu[%d] = %f\n", pts[i], pts[i], MinusK[pts[i]], pts[i], PlusK[pts[i]], pts[i], theta[pts[i]], pts[i], Nu[pts[i]]);
    }
   points[0] = 1;
   z[0]= n+1;
   printf("\nThe points on the Axis line are as follows:\n\n");
   printf("Points[0] = %d\n",points[0]);
   for(i=0;i<n-1;i++)
   {  points[i+1]= points[i] + z[i];
      printf("Points[%d] = %d\n",(i+1),points[i+1]);
      z[i+1]=z[i]-1;
   }
   MinusK[points[0]] = MinusK[pts[0]];
   PlusK[points[0]] = PlusK[pts[0]];
   theta[points[0]] = theta[pts[0]];
   Nu[points[0]] = Nu[pts[0]];
   printf("\n Flow properties of the points on the Axis line: \n");
   printf("\nMinusK[%d] = %f\tPlusK[%d] = %f\tTheta[%d] = %f\tNu[%d] = %f\n", points[0] ,
   MinusK[points[0]], points[0] , PlusK[points[0]], points[0] , theta[points[0]], points[0] ,
   Nu[points[0]]);
   for (i=1; i<n; i++)
   {  theta[points[i]]=0;
      Nu[points[i]]= Nu[pts[i]]+theta[pts[i]];
      MinusK[points[i]]= theta[points[i]] + Nu[points[i]];
      PlusK[points[i]]= theta[points[i]] - Nu[points[i]];
      printf("MinusK[%d]=%f\tPlusK[%d]=%f\tTheta[%d]=%f\tNu[%d]=%f\n",points[i],
      MinusK[points[i]],points[i],PlusK[points[i]],points[i],theta[points[i]],points[i],Nu[point
      s[i]]);
   }
   pnts[0]= n;
   y[0]=n;
   printf("\n\nThe points on the Last Right Running Characteristic line:\n\n"); printf("Pnts[0]= %d\n",pnts[0]);
   for(i=0;i<n-1;i++)
   {  pnts[i+1]= pnts[i] + y[i];
      printf("Pnts[%d] = %d\n",(i+1),pnts[i+1]);
      y[i+1]=y[i]-1;
   }
   printf("\n Flow properties of the points on the Last Right Running Characteristic line:\n");
   theta[pnts[0]] = t;
   Nu[pnts[0]] = t;
   MinusK[pnts[0]] = theta[pnts[0]] + Nu[pnts[0]];
   PlusK[pnts[0]] = theta[pnts[0]] - Nu[pnts[0]];
   printf("\nMinusK[ %d] = %f\tPlusK[ %d] = %f\tTheta[ %d] = %f\tNu[ %d] = %f\n",pnts[0],
   MinusK[pnts[0]],pnts[0], PlusK[pnts[0]],pnts[0], theta[pnts[0]],pnts[0], Nu[pnts[0]]);
   for (i=1; i<n;i++) 
   {  theta[pnts[i]]= dt*(n-k);
      Nu[pnts[i]] = - theta[points[i]] + Nu[points[i]] + theta[pnts[i]];
	  MinusK[pnts[i]]=theta[pnts[i]] + Nu[pnts[i]];
	  PlusK[pnts[i]]= theta[pnts[i]] - Nu[pnts[i]];
      printf("MinusK[%d] = %f\tPlusK[%d] = %f\tTheta[%d] = %f\tNu[%d] = %f\n",pnts[i],
      MinusK[pnts[i]],pnts[i], PlusK[pnts[i]], pnts[i],theta[pnts[i]], pnts[i],Nu[pnts[i]]);
      k++;
   }
   wallpts[0] = n+1;
   y[0] = n;
   printf("\n\nThe points on the Nozzle Wall:\n\n");
   printf("WallPts[0] = %d\n",wallpts[0]);
   for(i=0;i<n-1;i++)
   {  wallpts[i+1] = wallpts[i] + y[i];
      printf("WallPts[%d] = %d\n",(i+1),wallpts[i+1]);
      y[i+1] = y[i] - 1;
   }
   printf("\nFlow properties of the points on the Nozzle Wall: \n\n");
   for (i=0; i<n; i++)
   {  theta[wallpts[i]] = theta[pnts[i]];
      Nu[wallpts[i]] = Nu[pnts[i]];
      MinusK[wallpts[i]] = MinusK[pnts[i]];
      PlusK[wallpts[i]] = PlusK[pnts[i]];
      printf("MinusK[%d] = %f\tPlusK[%d] = %f\tTheta[%d] = %f\tNu[%d] = %f\n",wallpts[i],
      MinusK[wallpts[i]], wallpts[i], PlusK[wallpts[i]], wallpts[i],theta[wallpts[i]],
      wallpts[i],Nu[wallpts[i]]);
   }
   total = n+1;
   for(i=0; i<n-1; i++)
   total = total + (n-i);
   printf("\nTotal points = %d\n\n",total);
   for(i=1, p=1; i<=(total-n); i++, p++)
   {  for(j=0;j<n;j++)
      {  if(p == wallpts[j])
         {  jo = j;
            break;
         }
      }
      if(p == wallpts[j])
      {  intpts[i] = ++p;
         continue;
      }
      intpts[i] = p;
   }
   printf("\nThe Interior Points are:\n\n");
   for(i=0; i<(total-n); i++)
   {  intpts[i] = intpts[i+1];
      printf("Intpts[%d] = %d\n",i,intpts[i]);
   }
   printf("\n Flow properties of Interior points on Characteristic line 1 are:\n");
   for (i=0; i<n; i++)
   {  theta[intpts[i]] = theta[pts[i]];
      Nu[intpts[i]] = Nu[pts[i]];
      MinusK[intpts[i]] = MinusK[pts[i]];
      PlusK[intpts[i]] = PlusK[pts[i]];
      printf("\nMinusK[%d]=%f\tPlusK[%d]=%f\tTheta[%d]=%f\tNu[%d]=%f\n",intpts[i],
      MinusK[intpts[i]], intpts[i], PlusK[intpts[i]], intpts[i],theta[intpts[i]], intpts[i],Nu[intpts[i]]);
      for(i=0; i<n-1; i++)
      {  printf("\n\n Flow properties of Interior points on Characteristic line %d are:\n ",++q);
         theta[intpts[first]] = theta[points[i]];
         Nu[intpts[first]] = Nu [points[i]];
         MinusK[intpts[first]] = MinusK[points[i]];
         PlusK[intpts[first]] = PlusK[points[i]];
         first = first + (n-i);
         printf("\nMinusK[%d] = %f\tPlusK[%d] = %f\tTheta[%d] = %f\tNu[%d] = %f\n",intpts[first],
         MinusK[intpts[first]], intpts[first], PlusK[intpts[first]],
         intpts[first],theta[intpts[first]], intpts[first],Nu[intpts[first]]);
         second = second + (n-i);
         for(j=0;j<(n-x);j++)
         {  theta[intpts[second+j]] = (dt*(j+1));
            Nu[intpts[second+j]] = theta[intpts[second+j]] - PlusK[intpts[first]];
            MinusK[intpts[second+j]] = theta[intpts[second+j]] + Nu[intpts[second+j]]; PlusK
            [intpts[second+j]] = theta[intpts[second+j]] - Nu[intpts[second+j]];
            printf("\nMinusK[%d] = %f\tPlusK[%d] = %f\tTheta[%d] = %f\tNu[%d] = %f\n",intpts[second+j],MinusK[intpts[second+j]],intpts[second+j],PlusK[intpts[second+j]],
            intpts[second+j],theta[intpts[second+j]], intpts[second+j],Nu[intpts[second+j]]);
         }
         x=x+1;
      }
      printf("\n\n\n Flow properties of all Interior Points:\n");
      for(i=0; i<(total-n); i++)
      {  printf("\nMinusK[%d] = %f\tPlusK[%d] = %f\tTheta[%d] = %f\tNu[%d] = %f\n",intpts[i],
         MinusK[intpts[i]], intpts[i], PlusK[intpts[i]], intpts[i],theta[intpts[i]],
         intpts[i],Nu[intpts[i]]);
      }
      printf("\n\n\n Flow properties of the points on the Nozzle Wall:\n");
      for (i=0; i<n; i++)
      {  printf("\nMinusK[%d] = %f\tPlusK[%d] = %f\tTheta[%d] = %f\tNu[%d] = %f\n",wallpts[i],
         MinusK[wallpts[i]], wallpts[i], PlusK[wallpts[i]], wallpts[i],theta[wallpts[i]],
         wallpts[i],Nu[wallpts[i]]);
      }
   }
 }
