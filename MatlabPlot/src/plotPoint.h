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


#ifndef GRAPHPLOT_H
#define GRAPHPLOT_H

#include <vector>
#include <list>
#include <string>
#include <map>
#include <utility>
#include <initializer_list>
#include <fstream>

#include "common.h"
#include "linespec.h"

/*
 * 1. Markers are mostly the same (up to pt 13) except for terminal aqua.
 * 2. Only pngcairo and postscript export types support dashed lines.
 *    (jpeg, svg do not)
 * 3. Implement export for .png (pngcairo), .eps (postscript eps),
 *    .svg (no dashed), and monitor display (wxt or aqua)
 *
 * pngcairo: best
 * postscript: best
 * aqua: different markers no dashed
 * svg: no dashed
 */


namespace graph{

const std::string version = "0.1.0";

typedef std::vector<double> DataVector;

class PlotPoint
{
public:
    PlotPoint(unsigned mode=PDF);
    
    static PlotPoint* getInstance(unsigned mode=PDF);
    static PlotPoint* ints;
    static int instance;
    
    void xlabel(const std::string &label);
    void ylabel(const std::string &label);
    void title(const std::string &label);
    void legend(std::initializer_list<std::string> legendVec);
    void linespec(unsigned lineIndex, LineSpecInput lineSpec);
    void linespec(unsigned lineIndex, LineProperty property, std::string value);
    void linespec(unsigned lineIndex, LineProperty property, double value);
    void linespec(unsigned lineIndex,  std::string property, std::string value);
     
    void grid(bool flag);
    void plot(std::initializer_list<DataVector> il);
    void plot();
    
    void add(std::vector<double> &x);
    void print(const std::string &filenameExport);
    void exec(bool run_gnuplot=true);
    
    std::vector< std::vector<double> > vect1Plot ;
    std::string Marker;
    std::string MarkerValue;
    
    int size ;
    
    int vectCount; 
    int markerCount;

private:
    std::string filename;
    std::string labelX;
    std::string labelY;
    std::string labelTitle;
    std::vector<std::string>        legendVec;
    std::list<std::pair<unsigned, LineSpecInput>> lineSpecInput;
    std::vector<LineSpec>           lineSpec;
    std::list<std::string>          lineSpecAqua;
    std::list<std::string>          lineSpecCanvas;
    std::list<std::string>          lineSpecOther;
    unsigned nCurve;
    bool isGridded;
    //std::string filenameExport;

    bool flagScreen;
    bool flagHtml;
    bool flagPng;
    bool flagEps;
    bool flagPdf;
    bool flagSvg;

    bool existsAqua;
    bool existsWxt;
    bool existsCanvas;
    bool existsCairo;
    bool existsSvg;

    bool existsTerminal(std::string);
    void prepareLineSpec();

    //* plot curve .gp files
    void gpScreen(bool run_gnuplot);
    void gpPng(bool run_gnuplot);
    void gpEps(bool run_gnuplot);
    void gpPdf(bool run_gnuplot);
    void gpHtml(bool run_gnuplot);
    void gpSvg(bool run_gnuplot);
    void gpHeader(std::ofstream &fout);
    void gpCurve(std::ofstream &fout, const std::string &filename, bool run_gnuplot);
    
    
    std::string exportFileName();
    std::string prefixfilename();
};

}

#endif // GRAPHPLOT_H
