/*  (x)_n 
    vectorised on second argument
*/

void pochhammer(int *x, int n[], int *length, int answer[]){
  int i,j;
  
  for(j=0; j< *length; j++){
    answer[j]=1;
    for(i=0; i< n[j]; i++){
      answer[j] *= (*x)+i;
    }
  }
  return;
}
