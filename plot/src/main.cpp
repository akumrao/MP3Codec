#include "MatLabPlot.h"
#include "common.h"
using namespace std;
using namespace graph;

#ifndef M_PI
#define M_PI  3.14159265358979323846
#endif

double normal_pdf(double x, double mean, double sigma);







int main() {
    
    vector<double> x{1, 2, 3, 4, 5};
    vector<double> y {25, 15, 5, -5, -15};

    
    vector<double> x1 {1, 2, 3, 4, 5};
    vector<double> y2{1, 2, 3, 4, 8};
    
    try 
    {
        plot(x, y, LineStyle, "--" );
    
        plot(y);
    }
    catch(std::exception & e)
    {
      cout <<  e.what() << std::endl;
    }
    
    
    xlabel("x");
    ylabel("y"); // enhanced texts

    legend({"({/Symbol m}, {/Symbol s}) = (1, 1)", "({/Symbol m}, {/Symbol s}) = (2, 2)"});
    
    grid(true); // turn on grids

    
    try 
    {
        plot(x, y, "Color", "r", "LineWidth", 10, x, y2, "Color", "b", "LineWidth", 1);
    }
    catch(std::exception & e)
    {
      cout <<  e.what() << std::endl;
    }
    
   
    try
    {
    // Prepare sample data
    unsigned nPoint = 51;
    vector<double> x;
    x = linspace(0, 8, nPoint);
    vector<double> pdf11(nPoint);
    vector<double> pdf22(nPoint);
    double mu1 = 1;
    double mu2 = 2;
    double sigma1 = 1;
    double sigma2 = 2;

    for (unsigned i = 0; i < nPoint; ++i) {
        pdf11[i] = normal_pdf(x[i], mu1, sigma1);
        pdf22[i] = normal_pdf(x[i], mu2, sigma2);
    }


    xlabel("x");
    ylabel("pdf_{{/Symbol m},{/Symbol s}} (x)"); // enhanced texts

    legend({"({/Symbol m}, {/Symbol s}) = (1, 1)", "({/Symbol m}, {/Symbol s}) = (2, 2)"});
    
    grid(true); // turn on grids

    graph::plot(x, pdf11, "Color", "r", "LineWidth", 10, x, pdf22, "Color", "b", "LineWidth", 1);
        
    }
    catch(std::exception & e)
    {
      cout <<  e.what() << std::endl;
    }
}

double normal_pdf(double x, double mean, double sigma) {
    x -= mean;
    x /= sigma;
    return 1.0 / sqrt(2 * M_PI) / sigma * exp(-x * x / 2);
}

