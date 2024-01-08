/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Author:     Sami Byaruhanga 
Purpose:    implementing mandelbrot in OPENMP scheduling
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/


//INCLUDES + DEFINES *************************************************
#include <iostream>
#include <omp.h>
#include <vector>
#include <stdlib.h> 
#include <unistd.h>
#include<opencv4/opencv2/core/core.hpp>
#include<opencv4/opencv2/highgui/highgui.hpp>
// #include <opencv2/core/core.hpp> //for lab computers
// #include <opencv2/highgui/highgui.hpp>
// #include <opencv2/opencv.hpp>
using namespace std;
//************************************************* INCLUDES END HERE


//FUNCTIONS  ********************************************************
// struct to generate the mandelbrot set ----------------
typedef struct{
	double real;
	double imag;
	int colour;
}
complex_with_colour_t;
//-------------------------------------------------------

//-------------------------------------------------------
// processMandelbrotPixel 
// PURPOSE: To evaluate whether or not a pixel belongs to set
// INPUTS:  Max num of iterations to run recurrence, cutoff decision criterion, 
//          c i.e., point in complex plane
// RETURNS: c i.e., contains number of iterations in the recurrence
//------------------------------------------
void processMandelbrotPixel(int maxiter, double cutoff, complex_with_colour_t& c){
	complex_with_colour_t z; z.real = 0; z.imag = 0; //Setup initial recurrence variable z
	double tmp; 	//temp storage for the real part of z
	int count = 0; 	//number of iterations 
	double zabs; 	//|z|
	
	do 
	{//EVALUATE: z_{k+1} = z_{k}^2 + c;
		//z_{n+1} = z_{n}^2 + c;
		tmp = z.real*z.real - z.imag*z.imag;	//real part of z*z			
		z.imag = 2.0*z.real*z.imag;				//imag part of z*z
		z.real = tmp;
		
		z.real += c.real;						//add c			
		z.imag += c.imag;
		
		zabs = z.real*z.real + z.imag*z.imag;		//z*conj(z)
		count++;
		
	}
	while (zabs < cutoff && count < maxiter);
	c.colour = count;
}
//-------------------------------------------------------

//-------------------------------------------------------
// colorPixelFromScalar
// PURPOSE: Converts signle int colot into 3 ch img
// INPUT:   single double f [0,1]
// OUTPUT:  3 ch uchar coloyr can be used with opencv
//------------------------------------------
void colorPixelFromScalar(double f, cv::Vec3b& pixel){
	assert(f >= 0 - 1e-9 && f <= 1.0 + 1e-9);
	
	if (f < 0.03125) {pixel.val[2] = 59; pixel.val[1] = 76; pixel.val[0] = 192;}
	else if (f < 0.0625) {pixel.val[2] = 68; pixel.val[1] = 90; pixel.val[0] = 204;}
	else if (f < 0.09375) {pixel.val[2] = 77; pixel.val[1] = 104; pixel.val[0] = 215;}
	else if (f < 0.125) {pixel.val[2] = 87; pixel.val[1] = 117; pixel.val[0] = 225;}
	else if (f < 0.15625) {pixel.val[2] = 98; pixel.val[1] = 130; pixel.val[0] = 234;}
	else if (f < 0.1875) {pixel.val[2] = 108; pixel.val[1] = 142; pixel.val[0] = 241;}
	else if (f < 0.21875) {pixel.val[2] = 119; pixel.val[1] = 154; pixel.val[0] = 247;}
	else if (f < 0.25) {pixel.val[2] = 130; pixel.val[1] = 165; pixel.val[0] = 251;}
	else if (f < 0.28125) {pixel.val[2] = 141; pixel.val[1] = 176; pixel.val[0] = 254;}
	else if (f < 0.3125) {pixel.val[2] = 152; pixel.val[1] = 185; pixel.val[0] = 255;}
	else if (f < 0.34375) {pixel.val[2] = 163; pixel.val[1] = 194; pixel.val[0] = 255;}
	else if (f < 0.375) {pixel.val[2] = 174; pixel.val[1] = 201; pixel.val[0] = 253;}
	else if (f < 0.40625) {pixel.val[2] = 184; pixel.val[1] = 208; pixel.val[0] = 249;}
	else if (f < 0.4375) {pixel.val[2] = 194; pixel.val[1] = 213; pixel.val[0] = 244;}
	else if (f < 0.46875) {pixel.val[2] = 204; pixel.val[1] = 217; pixel.val[0] = 238;}
	else if (f < 0.5) {pixel.val[2] = 213; pixel.val[1] = 219; pixel.val[0] = 230;}
	else if (f < 0.53125) {pixel.val[2] = 221; pixel.val[1] = 221; pixel.val[0] = 221;}
	else if (f < 0.5625) {pixel.val[2] = 229; pixel.val[1] = 216; pixel.val[0] = 209;}
	else if (f < 0.59375) {pixel.val[2] = 236; pixel.val[1] = 211; pixel.val[0] = 197;}
	else if (f < 0.625) {pixel.val[2] = 241; pixel.val[1] = 204; pixel.val[0] = 185;}
	else if (f < 0.65625) {pixel.val[2] = 245; pixel.val[1] = 196; pixel.val[0] = 173;}
	else if (f < 0.6875) {pixel.val[2] = 247; pixel.val[1] = 187; pixel.val[0] = 160;}
	else if (f < 0.71875) {pixel.val[2] = 247; pixel.val[1] = 177; pixel.val[0] = 148;}
	else if (f < 0.75) {pixel.val[2] = 247; pixel.val[1] = 166; pixel.val[0] = 135;}
	else if (f < 0.78125) {pixel.val[2] = 244; pixel.val[1] = 154; pixel.val[0] = 123;}
	else if (f < 0.8125) {pixel.val[2] = 241; pixel.val[1] = 141; pixel.val[0] = 111;}
	else if (f < 0.84375) {pixel.val[2] = 236; pixel.val[1] = 127; pixel.val[0] = 99;}
	else if (f < 0.875) {pixel.val[2] = 229; pixel.val[1] = 112; pixel.val[0] = 88;}
	else if (f < 0.90625) {pixel.val[2] = 222; pixel.val[1] = 96; pixel.val[0] = 77;}
	else if (f < 0.9375) {pixel.val[2] = 213; pixel.val[1] = 80; pixel.val[0] = 66;}
	else if (f < 0.96875) {pixel.val[2] = 203; pixel.val[1] = 62; pixel.val[0] = 56;}
	else if (f < 1.0) {pixel.val[2] = 192; pixel.val[1] = 40; pixel.val[0] = 47;}
	else {pixel.val[2] = 180; pixel.val[1] = 4; pixel.val[0] = 38;}
}
//-------------------------------------------------------------
//************************************************ FUNCTIONS END HERE


//MAIN **************************************************************
int main(int argc, char**argv){
    //Initializations + CMD inputs + complex plane conversion
    if (argc < 9){
        std::cerr << "Command line parameters must include the following: \n"
            << "\t number of points in rows\n"
            << "\t number of points in cols\n"
            << "\t coordinate boundaries in real - xmin\n"  
            << "\t coordinate boundaries in real - xmax\n"  
            << "\t coordinate boundaries in imaginary - ymin\n" 
            << "\t coordinate boundaries in imaginary - ymax\n" 
            << "\t number of threads\n" 
            << "\t string type of scheduling to use either: dynamic, static, guided, runtime, auto" 
        << std::endl;
        return -1;
    }
    int ny = atoi(argv[1]);	//number of points in y (rows)
	int nx = atoi(argv[2]);	//number of points in x (cols)
	double xmin = strtod(argv[3], NULL);	//coorindate boundaries - x(real), y(img)
	double xmax = strtod(argv[4], NULL);
	double ymin = strtod(argv[5], NULL);
	double ymax = strtod(argv[6], NULL);
	double dx = (xmax - xmin)/(nx-1);	//Spacing -  recall that num or intervals is (num of points - 1)
	double dy = (ymax - ymin)/(ny-1);
    
    int nthreads = atoi(argv[7]);
    string scheduleName = argv[8];
    int chunkSize = 0;
    if(argv[9]!= NULL){chunkSize = atoi(argv[9]);}
    
    int maxiter = 1000; //Mandelbrot Set Parameters 
	double cutoff = 2;

    double serialTime = 0.0; //Timing 
    double parallelTime = 0.0;
    //--------------------------------------------------

    // Generate the points in the plane ----------------
    std::vector<std::vector<complex_with_colour_t> > points(ny);
    for (unsigned int iy = 0; iy < ny; iy++){ 
        points[iy].resize(nx);
        for (unsigned int ix = 0; ix < nx; ix++){
            points[iy][ix].real = xmin + ix*dx;
            points[iy][ix].imag = ymin + iy*dy;
            points[iy][ix].colour = 0;
        }
    }
    //----------------------------- END OF COMPLEX PLANE

    //Serial solution ---------------------------------
    double tstartSerial = omp_get_wtime();
    std::cout << "Serial code running ...\n" <<std::endl;
    //Create the set
    int serialCount = 0;
    for (unsigned int iy = 0; iy < ny; iy++){
        for (unsigned int ix = 0; ix < nx; ix++){
            processMandelbrotPixel(maxiter, cutoff, points[iy][ix]);
            serialCount += points[iy][ix].colour;
            
        }
    }
    //fill image
    cv::Mat mandelbrotSerial(ny, nx, CV_8UC3);
    for (unsigned int iy = 0; iy < ny; iy++){
        for (unsigned int ix = 0; ix < nx; ix++){
            double iteration_colouring = 1.0 - (double)points[iy][ix].colour/(double)maxiter;
            colorPixelFromScalar(iteration_colouring, mandelbrotSerial.at<cv::Vec3b>(iy,ix));				
        }
    }
    double tStopSerial = omp_get_wtime();
    serialTime = tStopSerial - tstartSerial;
    //Display the Serial image
    ostringstream converter;
    converter << "Serial Mandelbrot Set over x = " << xmin << ", " << xmax << " y = " << ymin << ", " << ymax <<".jpg";
    cv::imwrite(converter.str(),mandelbrotSerial);
    //--------------------------------------------------

	// Reset Image for Parallel ------------------------
    for (unsigned int iy = 0; iy < ny; iy++){
        for (unsigned int ix = 0; ix < nx; ix++){
            points[iy][ix].colour = 0;			
        }
    }
    //--------------------------------------------------
		
    //Parallel Solution -------------------------------
    omp_set_num_threads(nthreads);
    double tstartParallel = omp_get_wtime();

    //For different scheduling types ----------------------
    int parallelCount = 0;
    if(scheduleName == "dynamic"){//Can have chunksize
        if(chunkSize !=0){
            std::cout << "Dynamic scheduling with chunk size = " << chunkSize << " running ..." <<std::endl;
            #pragma omp parallel for schedule(dynamic, chunkSize)
            for (unsigned int iy = 0; iy < ny; iy++){
                int count = 0;
                for (unsigned int ix = 0; ix < nx; ix++){
                    processMandelbrotPixel(maxiter, cutoff, points[iy][ix]);
                    parallelCount += points[iy][ix].colour;
                }
            }
        }else{
            std::cout << "Dynamic scheduling running ..." <<std::endl;
            #pragma omp parallel for schedule(dynamic)
            for (unsigned int iy = 0; iy < ny; iy++){
                int count = 0;
                for (unsigned int ix = 0; ix < nx; ix++){
                    processMandelbrotPixel(maxiter, cutoff, points[iy][ix]);
                    parallelCount += points[iy][ix].colour;
                }
            }
        }
    }else if(scheduleName == "static"){//can have chunk size
        if(chunkSize !=0){
            std::cout << "static scheduling with chunk size = " << chunkSize << " running ..." <<std::endl;
            #pragma omp parallel for schedule(static, chunkSize)
            for (unsigned int iy = 0; iy < ny; iy++){
                int count = 0;
                for (unsigned int ix = 0; ix < nx; ix++){
                    processMandelbrotPixel(maxiter, cutoff, points[iy][ix]);
                    parallelCount += points[iy][ix].colour;
                }
            }
        }else{
            std::cout << "Static scheduling running ..." <<std::endl;
            #pragma omp parallel for schedule(static)
            for (unsigned int iy = 0; iy < ny; iy++){
                int count = 0;
                for (unsigned int ix = 0; ix < nx; ix++){
                    processMandelbrotPixel(maxiter, cutoff, points[iy][ix]);
                    parallelCount += points[iy][ix].colour;
                }
            }
        }
    }else if(scheduleName == "guided"){//can have chunk size
        if(chunkSize !=0){
            std::cout << "guided scheduling with chunk size = " << chunkSize << " running ..." <<std::endl;
            #pragma omp parallel for schedule(guided, chunkSize)
            for (unsigned int iy = 0; iy < ny; iy++){
                int count = 0;
                for (unsigned int ix = 0; ix < nx; ix++){
                    processMandelbrotPixel(maxiter, cutoff, points[iy][ix]);
                    parallelCount += points[iy][ix].colour;
                }
            }
        }else{
            std::cout << "Guided scheduling running ..." <<std::endl;
            #pragma omp parallel for schedule(guided)
            for (unsigned int iy = 0; iy < ny; iy++){
                int count = 0;
                for (unsigned int ix = 0; ix < nx; ix++){
                    processMandelbrotPixel(maxiter, cutoff, points[iy][ix]);
                    parallelCount += points[iy][ix].colour;
                }
            }
        }
    }else if(scheduleName == "runtime"){//No chunk size
        std::cout << "Runtime scheduling running ..." <<std::endl;
        #pragma omp parallel for schedule(runtime)
        for (unsigned int iy = 0; iy < ny; iy++){
            int count = 0;
            for (unsigned int ix = 0; ix < nx; ix++){
                processMandelbrotPixel(maxiter, cutoff, points[iy][ix]);
                parallelCount += points[iy][ix].colour;
            }
        }
    }else if(scheduleName == "auto"){//no chunk size
        std::cout << "Auto scheduling running ..." <<std::endl;
        #pragma omp parallel for schedule(auto)
        for (unsigned int iy = 0; iy < ny; iy++){
            int count = 0;
            for (unsigned int ix = 0; ix < nx; ix++){
                processMandelbrotPixel(maxiter, cutoff, points[iy][ix]);
                parallelCount += points[iy][ix].colour;
            }
        }
    }else{
        std::cerr << "You entered non-existent schedule name: " << scheduleName
            << " \nPlease use scheduling type either: dynamic, static, guided, runtime, auto" 
            << std::endl;
        return -1;
    }
    //--------------------------------------------------
    

    //Actually fill the image --------------------------
    cv::Mat mandelbrotParallel(ny, nx, CV_8UC3);
    for (unsigned int iy = 0; iy < ny; iy++){
        for (unsigned int ix = 0; ix < nx; ix++){
            double iteration_colouring = (double)points[iy][ix].colour/(double)maxiter; //NO 1-
            colorPixelFromScalar(iteration_colouring, mandelbrotParallel.at<cv::Vec3b>(iy,ix));				
        }
    }
    //--------------------------------------------------
    
    //Stop timer + print time + display results --------
    double tstopParallel = omp_get_wtime();
    parallelTime = tstopParallel - tstartParallel;
    ostringstream converter2;
    converter2 << scheduleName << " Mandelbrot Set over x = " << xmin << ", " << xmax << " y = " << ymin << ", " << ymax << ".jpg";
    cv::imwrite(converter2.str(), mandelbrotParallel);
    //--------------------------------------------------

    //Print Results ------------------------------------
    std::cout << "\nSerial Iterations =>" << serialCount << std::endl;
    std::cout << "Parallel iterations==> "<< parallelCount << std::endl;
        
    std::cout << "\nSerial Time = " << serialTime << " seconds." << std::endl;
    std::cout << "Parallel Time = " << parallelTime << " seconds." << std::endl;
    std::cout << "\nSpeedup = " << serialTime/parallelTime << std::endl;		
    //--------------------------------------------------
    
}
//**************************************************** MAIN ENDS HERE


/**STEPS TO SOLVE :::::::::::::::::::::::::::::::::::::::::::::::::::
 * 1. Use provided lecture MPI coded
 * 2. Get cmdline parameters
 * 3. do serial part and reset points 
 * 4. do OPENMP in parallel for and account for each schedule type
 * 5. print the speed up 
 * 
 * IMPROVEMENTS:
 * 1. Could improve code by writing the parallel part in helper function and call in main
 * 2. ...
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/