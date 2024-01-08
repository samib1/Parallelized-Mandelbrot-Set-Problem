# Parallelized Mandelbrot Set Problem

**Purpose:** To implement a shared memory version of the Mandelbrot set problem.

**Author:** Sami

**Program Files included:**  README.md, shared_mandelbrot.cpp

---


# Overview
The Mandelbrot Set problem is produced by applying a recurrence relation 𝑍𝑘+1=𝑍𝐾^2+𝐶 finite number of times. Then we color the point C based on the number of iterations it takes the recurrence relation to diverge. Divergence is quantified by checking |𝑍|^2 is greater than cutoff. 

Accelerated this problem by doing the processing of the Mandelbrot set in a parallel using various OpenMP scheduling types (namely; static, dynamic, guided, auto, runtime. Static, dynamic and guided can have chunk size parameters).

---

# How To Run The Program
## 1. Running 
1. compile with ``` make ```

2. Run command => ```./shared_mandelbrot rows cols xmin xmax ymin ymax nthreads scheduleName chunksize```.
Sample run command => ```./shared_mandelbrot 1024 1024 -1.5 0.5 -1 1 8 dynamic 5 >> ``` in cmdline.

3. Clean up with ``` make clean```

## 2. Sample output
1. Based off the parameters two images (shown below using the sample run command in 1 above) will be generated for the serial version and the parallel version of the mandelbrot set and the **speed up** will be printed on the terminal.

**Serial Output**

<img src = "https://github.com/samib1/Parallelized-Mandelbrot-Set-Problem/blob/main/Serial%20Mandelbrot%20Set%20over%20x%20%3D%20-1.5%2C%200.5%20y%20%3D%20-1%2C%201.jpg" width=50% height =50%> 

**Parallel Output**

<img src = "https://github.com/samib1/Parallelized-Mandelbrot-Set-Problem/blob/main/dynamic%20Mandelbrot%20Set%20over%20x%20%3D%20-1.5%2C%200.5%20y%20%3D%20-1%2C%201.jpg" width=50% height =50%>

---
## NOTES
1. To run this program you must c++, openCV, openMP installed.

___

**Have a wonderful day!**
