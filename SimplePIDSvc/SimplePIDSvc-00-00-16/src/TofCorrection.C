double tofsqrcov(double bg, vector<double> &data){
  

  if(bg<2)
    return data[0];
  else 
    return data[1];

  return data[0];
} 


double tofsigmascale(double x,  vector<double> &data){
  
  if(x<2)
    return data[2]*exp(data[3]*x)+data[4];
  else if(x<10) 
    return data[5];
  else
    return 1.0;
  
  return 1.0;
  
} 


double toflayer1scale(double x, double y,  vector<double> &data){
  
  double scale1=1.0;
  double scale1_cos=1.0;
  
  if(x<2){
    scale1= (data[6]*exp(data[7]*x)+data[8]*exp(data[9]*x)+data[10])/
      (data[11]*exp(data[12]*x)+data[13]*exp(data[14]*x)+data[15]);
    if(fabs(y)<=0.83)
      scale1_cos = (data[16]+data[17]*y+data[18]*y*y+data[19]*y*y*y)/(data[20]+data[21]*y+data[22]*y*y+data[23]*y*y*y);
  }
  else if(x<10){
    scale1= (data[24]*exp(data[25]*x)+data[26]*exp(data[27]*x)+data[28])/
      (data[29]*exp(data[30]*x)+data[31]*exp(data[32]*x)+data[33]);
    scale1_cos=1.0;
  }
  else{
    scale1=1.0;
    scale1_cos=1.0;
  }
  
    
  return scale1*scale1_cos;

}


double toflayer2scale(double x, double y,  vector<double> &data){
  
  double scale2=1.0;
  double scale2_cos=1.0;
  
  
  if(x<2){
    scale2= (data[34]*exp(data[35]*x)+data[36]*exp(data[37]*x)+data[38])/
      (data[39]*exp(data[40]*x)+data[41]*exp(data[42]*x)+data[43]);    
    if(fabs(y)<=0.83)
      scale2_cos = (data[44]+data[45]*y+data[46]*y*y+data[47]*y*y*y)/(data[48]+data[49]*y+data[50]*y*y+data[51]*y*y*y);
  }
  else if(x<10){
    scale2= (data[52]*exp(data[53]*x)+data[54]*exp(data[55]*x)+data[56])/
      (data[57]*exp(data[58]*x)+data[59]*exp(data[60]*x)+data[61]);    
    scale2_cos=1.0;
  }
  else{
    scale2=1.0;
    scale2_cos=1.0;
  } 
  
  return scale2*scale2_cos;
  
}



double mctofsqrcov(double bg,  vector<double> &mc){


  if(bg<2)
    return mc[0];
  else 
    return mc[1];
  
  return mc[0];
} 


double mctofsigmascale(double x, vector<double> &mc){
  
  if(x<2)
    return mc[2]*exp(mc[3]*x)+mc[4];
  else
    return mc[5];
  
  return 1.0;
  
} 


double mctoflayer1scale(double x, double y, vector<double> &mc){
  
  double scale=1.0;
  double scale_cos=1.0;
  if(x<2){
    scale= (mc[6]*exp(mc[7]*x)+mc[8]*exp(mc[9]*x)+mc[10])/
      (mc[11]*exp(mc[12]*x)+mc[13]*exp(mc[14]*x)+ mc[15]);
  }

  return scale*scale_cos;
}


double mctoflayer2scale(double x, double y, vector<double> &mc){
  
  double scale=1.0;
  double scale_cos=1.0;

  if(x<2){
    scale= (mc[16]*exp(mc[17]*x)+mc[18]*exp(mc[19]*x)+mc[20])/
      (mc[21]*exp(mc[22]*x)+mc[23]*exp(mc[24]*x)+ mc[25]);
  }

  return scale*scale_cos;
}
