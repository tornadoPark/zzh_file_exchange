
  temp1=16.3*faceArea(18,10,1)* ( T(17,10)-T(18,10))/dx;
  temp2=16.3* faceArea(18,10,2)*( T(19,10)-T(18,10))/dx;
  temp3=16.3* faceArea(18,10,3)*( T(18,9)-T(18,10))/dr;
  temp4= 14.7023* faceArea(18,10,4)*( T(18,11)-T(18,10))/dr;
  temp5= volumetotal(18,10)*(T(18,10)-TLast(18,10))*8000*1260/dt;
  temp1+temp2+temp3+temp4