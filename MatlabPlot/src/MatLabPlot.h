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
#include "plotPoint.h"

using namespace std;

namespace graph{

void xlabel( std::string l);

void ylabel( std::string l);

void legend( std::initializer_list<std::string> l );

void grid( bool l);


class mytype {
public:
     mytype(int x) {
        vecType = false;
        str =to_string(x);
    }
     
    mytype(double x) {
       vecType = false;
       str =to_string(x);

    }
    mytype(string x) {
        vecType = false;
        str =x;
    }

    mytype(char *x) {
        vecType = false;
        cout << "char *x " << x;
    }
     void fill (auto &x) {
        vecType = true;
        for (auto i : x)
           vects.push_back( i);
    }

    mytype(list<double> &x) {
        fill(x);
    }
    
    mytype(list<float> &x) {
        fill(x);
    }
     
    mytype(list<int> &x) {
        fill(x);
    }
    mytype(vector<double> &x) {
        fill(x);
    }
     mytype(vector<int> &x) {
        fill(x);
    }

     mytype(vector<float> &x) {
        fill(x);
    }
    vector<double> vects;
    string str;
    bool vecType;
    
};

template<typename T>
void plotgraph(T t) {
    
    PlotPoint *curvePlot = PlotPoint::getInstance();
    
   // curvePlot->size = size;

    int &vectCount = curvePlot->vectCount;
    int &markerCount = curvePlot->markerCount;

    string &Marker = curvePlot->Marker;
    string &MarkerValue = curvePlot->MarkerValue;
    
    mytype d(t);
    if (d.vecType) {
        if (vectCount % 2 == 0) {
            curvePlot->add( d.vects);
        }
        if (vectCount % 2 == 1) {
            curvePlot->add( d.vects);
        }

        ++vectCount;
    } else {
        if (vectCount < 2) {
            cout << " First & Second paramter should be vector or list";
            return;
        } else {
            if (markerCount % 2 == 0) {
                Marker = d.str;
            }
            if (markerCount % 2 == 1) {
                MarkerValue = d.str;

               // curvePlot->plot({vect1, vect2});
                if( vectCount % 2 != 0)
                {
                    throw " some internal error contact developer";
                }
                    
                curvePlot->linespec(vectCount/2, Marker, MarkerValue );
            }

            ++markerCount;

        }
    }

}


template<typename T>
void plot(T t) {
    
    PlotPoint *curvePlot = PlotPoint::getInstance();
    
    int &vectCount = curvePlot->vectCount;
    int &markerCount = curvePlot->markerCount;
    
   
    if(vectCount || markerCount )
    plotgraph( t);
    
    else
    {
        mytype d(t);

        if( d.vecType)
        {

             vector<double> x;
             for( int i=0 ; i < d.vects.size() ; ++i)
                 x.push_back(i);

            curvePlot->add(  x);
            curvePlot->add(  d.vects);

        }
        else{
          cout << " First paramter should be vector or list"  ; 

        }
    }
    
    curvePlot->plot();
    
    curvePlot->exec();
    
    vectCount = 0;
    markerCount =0;
    
    curvePlot->vect1Plot.clear();
    
}

template<typename T, typename... R>
void plot(T t, R... rest) {
    
   //int size = sizeof...(rest);
    plotgraph( t); 
    plot(rest...);
}

}//End Namespace

#endif /* MATLABPLOT_H */

