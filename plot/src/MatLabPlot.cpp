/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   MatLabPlot.cpp
 * Author: Arvind Umrao
 * Email akumrao@yahoo.com  akumrao@gmail.com
 * Created on May 8, 2018, 12:00 AM
 */

#include "MatLabPlot.h"


using namespace std;

namespace graph{
    

void xlabel( string l)
{   PlotPoint *curvePlot = PlotPoint::getInstance();
    curvePlot->xlabel(l); 
}
void ylabel( string l)
{   
    PlotPoint *curvePlot = PlotPoint::getInstance();
    curvePlot->ylabel(l); 
} 

void legend( std::initializer_list<string> l )
{   PlotPoint *curvePlot = PlotPoint::getInstance();
    curvePlot->legend(l); 
} 
void grid( bool l)
{   PlotPoint *curvePlot = PlotPoint::getInstance();
    curvePlot->grid(l); 
} 



}//End Namespace