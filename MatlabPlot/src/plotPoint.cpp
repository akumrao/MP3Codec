#include "plotPoint.h"

#include<algorithm>
#include<fstream>
#include<stdexcept>
#include<cstdlib>
#include<iostream>
#include<sstream>
#include<cmath>

using namespace std;

namespace graph {

  std::map<std::string,LineProperty>LineSpecmap = { {"LineStyle", LineStyle}, {"LineWidth",LineWidth}, {"Marker", Marker}, {"MarkerSize", MarkerSize}, {"Color", Color} };

    
    int PlotPoint::instance = 0;

  PlotPoint* PlotPoint::ints = nullptr;
  PlotPoint*  PlotPoint::getInstance(unsigned mode)
  {
      if( ints == NULL )
      if( ints == NULL )
      {
          // make it thread safe in future if more call to plot
          ints = new PlotPoint(mode);
       }
     // ++instance;
      return ints;
  }
  
    PlotPoint::PlotPoint(unsigned mode)
    : filename("plot"),
    labelX(),
    labelY(),
    labelTitle(),
    legendVec(0),
    lineSpecInput(),
    lineSpec(),
    lineSpecAqua(),
    lineSpecCanvas(),
    lineSpecOther(),
    nCurve(0),
    isGridded(false)

    {
        vectCount = 0;
        markerCount =0;

        vect1Plot.empty();
        
         //* Set up modes
        flagScreen = (SCREEN & mode) ? true : false;
        flagPng = (PNG & mode) ? true : false;
        flagEps = (EPS & mode) ? true : false;
        flagPdf = (PDF & mode) ? true : false;
        flagHtml = (HTML & mode) ? true : false;
        flagSvg = (SVG & mode) ? true : false;
        
        //* Test if terminal exists
        if(flagScreen)
        this->existsAqua   = this->existsTerminal("aqua");
        //this->existsWxt    = this->existsTerminal("wxt");
        
        this->existsCairo = this->existsTerminal("cairo");
        //this->existsCanvas = this->existsTerminal("canvas");
        //
        if(flagSvg)
        this->existsSvg    = this->existsTerminal("svg");
     
    }

    void PlotPoint::xlabel(const std::string &label) {
        this->labelX = label;
    }

    void PlotPoint::ylabel(const std::string &label) {
        this->labelY = label;
    }

    void PlotPoint::title(const std::string &label) {
        this->labelTitle = label;
    }

    void PlotPoint::legend(std::initializer_list<string> legendVec) {
        this->legendVec.resize(legendVec.size());
        std::copy(legendVec.begin(), legendVec.end(), this->legendVec.begin());
    }

    void PlotPoint::linespec(unsigned lineIndex, LineSpecInput lineSpec) {
        //* Store but no process

        if (lineIndex <= 0) {
            throw out_of_range("Line index must be a positive integer");
        }

        this->lineSpecInput.push_back({lineIndex, lineSpec});
    }

   void PlotPoint::linespec(unsigned lineIndex,  string property, string value) {
       
       
       linespec( lineIndex, LineSpecmap[property] ,  value);
    }
        
    void PlotPoint::linespec(unsigned lineIndex, LineProperty property, string value) {
        linespec(lineIndex,{
            {property, value}});
    }

    void PlotPoint::linespec(unsigned lineIndex, LineProperty property, double value) {
        linespec(lineIndex,{
            {property, to_string(value)}});
    }

    void PlotPoint::grid(bool flag) {
        this->isGridded = flag;
    }

    void PlotPoint::add( vector<double> &x) {
        vect1Plot.push_back(x);
    }
    
    void PlotPoint::plot()
    {
        if (vect1Plot.size() % 2) {
            throw length_error("Arguements must be even number of data vectors");
        }

        //* check if columns are of equal lengths
        ofstream fout((prefixfilename() + string(".dat")).c_str());
        this->nCurve = 0;
        for (auto it = vect1Plot.begin(); it != vect1Plot.end(); ++it) {
            auto itEven = it++;
            if (it->size() != itEven->size()) {
                throw length_error("Pairwise data vectors must have the same lengths");
            }
            const DataVector &x = *itEven;
            const DataVector &y = *it;

            fout << "# Curve " << (this->nCurve)++ << endl;
            for (unsigned i = 0; i < x.size(); ++i) {
                fout << x[i] << "," << y[i] << endl;
            }
            fout << endl << endl;
        }
        fout.close();
        
        vect1Plot.clear();
    }
    
    void PlotPoint::plot(initializer_list<DataVector> il) {
        //* Take Matlab-like commands but only store data
        //* Actually plotting happens at PlotPoint::show()
 
       
        if (il.size() % 2) {
            throw length_error("Arguements must be even number of data vectors");
        }

        //* check if columns are of equal lengths
        ofstream fout((prefixfilename() + string(".dat")).c_str());
        this->nCurve = 0;
        for (auto it = il.begin(); it != il.end(); ++it) {
            auto itEven = it++;
            if (it->size() != itEven->size()) {
                throw length_error("Pairwise data vectors must have the same lengths");
            }
            const DataVector &x = *itEven;
            const DataVector &y = *it;

            fout << "# Curve " << (this->nCurve)++ << endl;
            for (unsigned i = 0; i < x.size(); ++i) {
                fout << x[i] << "," << y[i] << endl;
            }
            fout << endl << endl;
        }
        fout.close();
    }

    void PlotPoint::print(const string &filename) {
        this->filename = filename;
    }

    string PlotPoint::exportFileName() {
        return this->filename+ "export" + std::to_string(instance);
    }

    string PlotPoint::prefixfilename() {
        return this->filename + std::to_string(instance);
    }

        
    void PlotPoint::exec(bool run_gnuplot) {
        //* Check if there are data
        if (this->nCurve == 0) {
            return;
        }            //* Check if legend size is zero
        else if (this->legendVec.size() == 0) {
            this->legendVec.resize(this->nCurve);
            for (unsigned i = 0; i<this->nCurve; ++i) {
                this->legendVec[i] = to_string(i + 1);
            }
        }            //* Check if legend size matches data pair number. If not, pad it.
        else if (this->legendVec.size() != this->nCurve) {
            //throw length_error("Legends must match the number of data vector pairs");
            unsigned nLegend = this->legendVec.size();
            this->legendVec.resize(this->nCurve);
            for (unsigned i = nLegend; i < nCurve; ++i) {
                this->legendVec[i] = "Data " + to_string(i + 1);
            }
        }

        prepareLineSpec();

        if (this->flagScreen) {
            gpScreen(run_gnuplot);
        }
        if (this->flagPng) {
            gpPng(run_gnuplot);
        }
        if (this->flagEps) {
            gpEps(run_gnuplot);
        }
        if (this->flagPdf) {
            gpPdf(run_gnuplot);
        }
        if (this->flagHtml) {
            gpHtml(run_gnuplot);
        }
        if (this->flagSvg) {
            gpSvg(run_gnuplot);
        }
        ++instance;
        vectCount = 0;
        markerCount =0;

        vect1Plot.empty();
        
    }

    bool PlotPoint::existsTerminal(string terminalName) {
        bool result = false;

        string commandTest = "gnuplot -e \"set print '" + prefixfilename()
                + "-exists-" + terminalName + "'; if (strstrt(GPVAL_TERMINALS, '" + terminalName + "')) print 1; else print 0\"";
        std::system(commandTest.c_str());

        ifstream fin((prefixfilename() + "-exists-" + terminalName).c_str());
        fin >> result;
        fin.close();

        return result;
    }

    void PlotPoint::prepareLineSpec() {
        LineSpec::resetLineCount();
        this->lineSpec.resize(nCurve);

        for (auto it = lineSpecInput.begin(); it != lineSpecInput.end(); ++it) {
            unsigned lineIndex = it->first;
            if (lineIndex > nCurve) {
                continue;
            }
            for (auto itProperty = it->second.begin(); itProperty != it->second.end(); ++itProperty) {
                this->lineSpec[lineIndex - 1].set(*itProperty);
            }
        }
    }

    void foutGridSetting(ofstream &fout, TerminalType tt) {
        fout << "set grid lc rgb '" << LineSpec::gridColor << "' lw 1 lt " << LineSpec::getGridLineType(tt) << endl;
    }

    void PlotPoint::gpScreen(bool run_gnuplot) {
        //* Generate gnuplot batch file
        string filename = prefixfilename() + ".gp";
        ofstream fout(filename.c_str());

        gpHeader(fout);

        //* Set terminal and line style
        if (existsAqua) {
            fout << "set terminal aqua dashed enhanced" << endl;
            if (this->isGridded) {
                foutGridSetting(fout, TERM_AQUA);
            }
            for (unsigned i = 0; i<this->lineSpec.size(); ++i) {
                fout << this->lineSpec[i].toStringAqua() << endl;
            }
        } else if (existsWxt) {
            fout << "set terminal wxt dashed enhanced" << endl;
            if (this->isGridded) {
                foutGridSetting(fout, TERM_WXT);
            }
            for (unsigned i = 0; i<this->lineSpec.size(); ++i) {
                fout << this->lineSpec[i].toStringWxtCairoSvg() << endl;
            }
        } else {
            fout << "# No supported display terminal found. Line styles may be not accurate." << endl;
            if (this->isGridded) {
                foutGridSetting(fout, TERM_OTHER);
            }
            for (unsigned i = 0; i<this->lineSpec.size(); ++i) {
                fout << this->lineSpec[i].toStringWxtCairoSvg() << endl;
            }
        }
        gpCurve(fout, filename, run_gnuplot);

        fout.close();

    }

    void PlotPoint::gpPng(bool run_gnuplot) {
        //* Generate gnuplot batch file
        string filename = prefixfilename() + "-png.gp";
        string filenameExport = exportFileName() + ".png";
        ofstream fout(filename.c_str());

        gpHeader(fout);

        //* Set terminal and line style
        if (existsCairo) {
            fout << "set terminal pngcairo dashed enhanced" << endl;
            fout << "set output '" << filenameExport << "'" << endl;
            if (this->isGridded) {
                foutGridSetting(fout, TERM_CAIRO);
            }
            for (unsigned i = 0; i<this->lineSpec.size(); ++i) {
                fout << this->lineSpec[i].toStringWxtCairoSvg() << endl;
            }
        } else {
            fout << "# Cairo terminal not found. Default png terminal used instead." << endl
                    << "# Line styles may be not accurate." << endl;
            fout << "set terminal png dashed enhanced" << endl;
            fout << "set output '" << filenameExport << "'" << endl;
            if (this->isGridded) {
                foutGridSetting(fout, TERM_OTHER);
            }
            for (unsigned i = 0; i<this->lineSpec.size(); ++i) {
                fout << this->lineSpec[i].toStringWxtCairoSvg() << endl;
            }
        }
        gpCurve(fout, filename, run_gnuplot);

        fout.close();
    }

    void PlotPoint::gpEps(bool run_gnuplot) {
        //* Generate gnuplot batch file
        string filename = prefixfilename() + "-eps.gp";
        string filenameExport = exportFileName() + ".eps";
        ofstream fout(filename.c_str());

        gpHeader(fout);

        //* Set terminal and line style
        string gridColor = "#cccccc";

        if (existsCairo) {
            fout << "set terminal epscairo transparent color dashed enhanced" << endl;
            fout << "set output '" << filenameExport << "'" << endl;
            if (this->isGridded) {
                foutGridSetting(fout, TERM_CAIRO);
            }
            for (unsigned i = 0; i<this->lineSpec.size(); ++i) {
                fout << this->lineSpec[i].toStringWxtCairoSvg() << endl;
            }
        } else {
            fout << "# Cairo terminal not found. Postscript terminal used instead." << endl
                    << "# Line styles may be not accurate." << endl;
            fout << "set terminal postscript eps color colortext dashed" << endl;
            fout << "set output '" << filenameExport << "'" << endl;
            if (this->isGridded) {
                foutGridSetting(fout, TERM_OTHER);
            }
            for (unsigned i = 0; i<this->lineSpec.size(); ++i) {
                fout << this->lineSpec[i].toStringWxtCairoSvg() << endl;
            }
        }

        gpCurve(fout, filename, run_gnuplot);

        fout.close();
    }

    void PlotPoint::gpPdf(bool run_gnuplot) {
        //* Generate gnuplot batch file
        string filename = prefixfilename() + "-pdf.gp";
        string filenameExport = exportFileName() + ".pdf";
        ofstream fout(filename.c_str());

        gpHeader(fout);

        //* Set terminal and line style
        if (existsCairo) {
            fout << "set terminal pdfcairo transparent color dashed enhanced" << endl;
            fout << "set output '" << filenameExport << "'" << endl;
            if (this->isGridded) {
                foutGridSetting(fout, TERM_CAIRO);
            }
            for (unsigned i = 0; i<this->lineSpec.size(); ++i) {
                fout << this->lineSpec[i].toStringWxtCairoSvg() << endl;
            }
        }

        gpCurve(fout, filename, run_gnuplot);

        fout.close();
    }

    void PlotPoint::gpHtml(bool run_gnuplot) {
        //* Generate gnuplot batch file
        string filename = prefixfilename() + "-html.gp";
        string filenameExport = exportFileName() + ".html";
        ofstream fout(filename.c_str());

        gpHeader(fout);

        //* Set terminal and line style
        if (existsCairo) {
            fout << "set terminal canvas dashed enhanced" << endl;
            fout << "set output '" << filenameExport << "'" << endl;
            if (this->isGridded) {
                foutGridSetting(fout, TERM_CANVAS);
            }
            for (unsigned i = 0; i<this->lineSpec.size(); ++i) {
                fout << this->lineSpec[i].toStringWxtCairoSvg() << endl;
            }
        }

        gpCurve(fout, filename, run_gnuplot);

        fout.close();
    }

    void PlotPoint::gpSvg(bool run_gnuplot) {
        //* Generate gnuplot batch file
        string filename = prefixfilename() + "-svg.gp";
        string filenameExport = exportFileName() + ".svg";
        ofstream fout(filename.c_str());

        gpHeader(fout);

        //* Set terminal and line style
        if (existsCairo) {
            fout << "set terminal svg dashed enhanced" << endl;
            fout << "set output '" << filenameExport << "'" << endl;
            if (this->isGridded) {
                foutGridSetting(fout, TERM_SVG);
            }
            for (unsigned i = 0; i<this->lineSpec.size(); ++i) {
                fout << this->lineSpec[i].toStringWxtCairoSvg() << endl;
            }
        }

        gpCurve(fout, filename, run_gnuplot);

        fout.close();

    }

    void PlotPoint::gpHeader(ofstream &fout) {
        fout << "# Gnuplot script file" << endl;
        fout << "# Automatically generated by plot Ver. " << version << endl;
        fout << "set datafile separator ','" << endl;
    }

    void PlotPoint::gpCurve(ofstream &fout, const string &filename, bool run_gnuplot) {
        fout << "set style increment userstyle" << endl;
        fout << "set autoscale" << endl;
        fout << "unset log" << endl;
        fout << "unset label" << endl;
        fout << "set xtic auto" << endl;
        fout << "set ytic auto" << endl;
        fout << "set title \"" << this->labelTitle << "\"" << endl;
        fout << "set xlabel \"" << this->labelX << "\"" << endl;
        fout << "set ylabel \"" << this->labelY << "\"" << endl;
        fout << "plot ";

        for (unsigned i = 0; i<this->nCurve; ++i) {

            string filename = prefixfilename() + ".dat";
            fout << "'" << filename << "' index " << i
                    << " title '" << this->legendVec[i]
                    << "' with ";

            if (this->lineSpec[i].isPointOnly()) {
                fout << "points";
            } else {
                fout << "linespoints";
            }
            fout << ", ";
        }
        fout << endl;

        if (run_gnuplot) {
            std::system(("gnuplot " + filename).c_str());
        }

    }



}
