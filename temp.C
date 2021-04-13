#include "DT_Classifier.h"
#include <stdio.h>
#include<fstream>
#include <vector>
#include<iostream>

int main()
{
const int32_t features_length = 6;

 const float features[features_length] = {300, 2000, 500, 0.125, 3.01,3.23};//{4.46000000e+02, 2.54778044e+03, 1.50210045e+02, 5.89572174e-02, 3.45358872e+00, 3.53086019e+00};

const int32_t predicted_class = DT_Classifier_predict(features, features_length);
 std::cout<< predicted_class<<std::endl; // 0= pions, 1 = electrons
 return 0;
}
