/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MatLabPlot.h
 * Author: root
 *
 * Created on May 8, 2018, 12:00 AM
 */

#ifndef MATLABPLOT_H
#define MATLABPLOT_H
#include <iostream>
#include <string>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>


namespace graph{

void xlabel( std::string l);

void ylabel( std::string l);

void legend( std::initializer_list<std::string> l );

void grid( bool l);



template<typename T>
void plotgraph(T t) ;
    
 


template<typename T>
void plot(T t);



template<typename T, typename... R>
void plot(T t, R... rest) ;

}//End Namespace

#endif /* MATLABPLOT_H */

