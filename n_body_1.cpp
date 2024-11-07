// n_body_1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <array>
#include <vector>
#include <string_view>

#include <iostream>
#include <fstream>
#include <cmath>

#include "progressbar.hpp"

// constant variables
constexpr double g{ 6.67430e-11 }; // value of gravitational constant (m^3 s^-2 kg^-1)
//constexpr double au{ 149597870700 }; // value of an astronomical unit (m)
constexpr double pi{ 3.141592653589793 }; // value of pi
constexpr int bodies{ 3 }; // number of interacting bodies in system

// definition of constant arrays
constexpr std::array<std::string_view, bodies> name{ "Sun", "Earth", "Moon" };//, "Jupiter", "Saturn", "Uranus", "Neptune" }; // string_view names of each body in system
constexpr std::array<double, bodies> gMass{ 1.32712440041279419e20, 398600435507000, 4902800118000 };// , 126712764100000000, 37940584841800000, 5794556400000000, 6836527100580000 }; // masses of each body in system multiplied by the graviational constant (m^3 s^-2)
std::array<double, bodies> mass{ gMass[0] / g, gMass[1] / g, gMass[2] / g };//, gMass[3] / g, gMass[4] / g, gMass[5] / g, gMass[6] / g }; // masses of each body in system (kg)
std::array<int, bodies> parent{ 0, 0, 1 };//, 0, 0, 0, 0 }; // the index of the parent of this body

// functions for analytics
// convert degrees input to radians
static double degToRad(double deg)
{
    return deg * (pi / 180);
}

// takes two 1D vector arrays and sums each element to each other
static std::vector<double> parentHandle(std::vector<double> (&parent), std::vector<double> (&child))
{
    if (parent.size() != child.size())
    {
        return { 0, 0, 0 };
    }

    std::vector<double> temp{ 0, 0, 0 };

    for (int i = 0; i < parent.size(); ++i)
    {
        temp[i] = parent[i] + child[i];
    }

    return temp;
}

static double eucDist(std::vector<std::vector<double>> (&d), int i, int j) // calculates the euclidian distance between two bodies
{
    return hypot(hypot(d[i][0]- d[j][0], d[i][1] - d[j][1]), d[i][2] - d[j][2]);
}

static double lengthVect(std::vector<std::vector<double>> (&d), int i)
{
    return sqrt( pow(d[i][0],2) + pow(d[i][1], 2) + pow(d[i][2], 2) );
}

// calculates the total angular momentum of the system
static double angMome(std::vector<std::vector<double>> (&d), std::vector<std::vector<double>> (&v), std::array<double, bodies> (&m), int n)
{
    double count = 0;

    for (int i = 0; i < n; ++i)
    {
        count += sqrt( pow(lengthVect(d, i), 2) * pow(lengthVect(v, i) * m[i], 2) - pow( d[i][0]*v[i][0]*m[i] + d[i][1] * v[i][1] * m[i] + d[i][2] * v[i][2] * m[i] , 2));
    }

    return count;
}

// calculates the total energy of a body relative to its parent body
static double totalEner(std::vector<std::vector<double>>(&d), std::vector<std::vector<double>>(&v), std::array<double, bodies>(&m), int i, int j)
{
    if (i == j)
    {
        return 0;
    }

    double totV{ 0 };
    double totR{ 0 };

    for (int k = 0; k < 3; ++k)
    {
        totV += pow(v[i][k] - v[j][k], 2);
        totR += pow(d[i][k] - d[j][k], 2);
    }

    return m[i] * m[j] * (0.5 * totV / (m[i] + m[j]) - g/sqrt(totR));
}

// a boolean function that takes in the total energy of a body and returns its bound status to its parent
static bool bound(double totE)
{
    if (totE < 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// calculates the angular momentum between a body and its parent
static double totalMome(std::vector<std::vector<double>>(&d), std::vector<std::vector<double>>(&v), std::array<double, bodies>(&m), int i, int j)
{
    double totL = 0;

    double ang1{ (d[i][1] - d[j][1]) * (v[i][2] - v[j][2]) - (d[i][2] - d[j][2]) * (v[i][1] - v[j][1]) };
    double ang2{ (d[i][2] - d[j][2]) * (v[i][0] - v[j][0]) - (d[i][1] - d[j][1]) * (v[i][2] - v[j][2]) };
    double ang3{ (d[i][0] - d[j][0]) * (v[i][1] - v[j][1]) - (d[i][1] - d[j][1]) * (v[i][0] - v[j][0]) };

    totL = sqrt(pow(ang1, 2) + pow(ang2, 2) + pow(ang3, 2)) * m[i] * m[j] / (m[i] + m[j]);

    return totL;
}

// calculates the eccentricity of a body around its parent
static double eccentricity(std::array<double, bodies>(&m), double totL, double totE, int i, int j) 
{
    double ecc{ 2 * (m[i] + m[j]) * pow(totL, 2) * totE / (pow(g, 2) * pow(m[i] * m[j], 3))};

    return sqrt(1 + ecc);
}

// functions for orbital mechanics
// function takes in mass (of parent), eccentricity, semi-major axis, argument of periapsis, inclination, longitude of ascending node, and true anomaly and returns a planets cartesian position and velocity at that point
static std::vector<std::vector<double>> kepToCar(double m, double e, double a, double ap, double i, double lan, double v)
{
    // converting degree inputs to radians
    ap = degToRad(ap);
    i = degToRad(i);
    lan = degToRad(lan);
    v = degToRad(v);

    double p{ a * (1 - pow(e,2)) };

    double r{ p / (1 + e * cos(v)) }; // calculates radius of orbit
    double h{ sqrt( g * m * p ) }; // calculated specific angular momentum of orbiting body

    double cLan{ cos(lan) };
    double sLan{ sin(lan) };
    double cApV{ cos(ap + v) };
    double sApV{ sin(ap + v) };
    double cI{ cos(i) };
    double sI{ sin(i) };

    double x{ r * (cLan * cApV - sLan * sApV * cI) };
    double y{ r * (sLan * cApV + cLan * sApV * cI) };
    double z{ r * sI * sApV };

    double xDot{ x * h * e * sin(v) / (r * p) - (h / r) * (cLan * sApV + sLan * cApV * cI) };
    double yDot{ y * h * e * sin(v) / (r * p) - (h / r) * (sLan * sApV - cLan * cApV * cI) };
    double zDot{ z * h * e * sin(v) / (r * p) + (h / r) * (sI * cApV) };

    return { { x, y, z }, { xDot, yDot, zDot } };
}

// functions for advancing the simulation
static void updateDisp(std::vector<std::vector<double>> (&d), std::vector<std::vector<double>> (&v), double step, int n) // updates the location of each interacting body 
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            d[i][j] = d[i][j] + v[i][j] * step / 2;
        }
    }
}

static void updateVelo(std::vector<std::vector<double>> (&v), std::vector<std::vector<double>> (&a), double step, int n) // updates the velocity of each interacting body 
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            v[i][j] = v[i][j] + a[i][j] * step;
        }
    }
}

static void updateAcce(std::vector<std::vector<double>> (&d), std::vector<std::vector<double>> (&a), std::array<double, bodies> (&m), int n) // updates the acceleration of each interacting body 
{
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            a[i][j] = 0;
        }
    }
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            if (j != i)
            {
                for (int k = 0; k < 3; ++k)
                {
                    a[i][k] -= (g * m[j] * (d.at(i).at(k) -
                        d[j][k])) / pow(eucDist(d, i, j), 3);
                }
            }
        }
    }
}

int main()
{

    // definition of motion arrays
    std::vector<std::vector<double>> disp{ {{0, 0, 0}, {149.6e9, 0, 0}, {149.984e9, 0, 0}, {778.5e9, 0, 0}, {1432e9, 0, 0}, {2867e9, 0, 0}, {4515e9, 0, 0}} }; // initial displacement of each body (m)
    std::vector<std::vector<double>> velo{ {{0, 0, 0}, {0, 29.8e3, 0}, {0, 30.8e3, 0}, {0, 13.1e3, 0}, {0, 9.7e3, 0}, {0, 6.8e3, 0}, {0, 5.4e3, 0}} }; // initial velocity of each body (m s^-1)
    std::vector<std::vector<double>> acce{ {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}, {0, 0, 0}} }; // initial acceleration of each body (m s^-2)

    // time variables
    double time{ 0 }; // total time of simulation (s)
    double timestep{ 100 }; // timestep of each simulation calculation (s)

    // numerical variables
    //int subBodies {}; // number of non-interacting bodies in system

    // changing the properties of earth to match given parameters
    std::vector<std::vector<double>> earthNew{ kepToCar(mass[0], 0.017, 149.6e9, 115.54, 0, 0, 0) };
    disp[1] = earthNew[0];
    velo[1] = earthNew[1];

    std::vector<std::vector<double>> moonNew{ kepToCar(mass[1], 0.055, 0.384e9, 0, 5.1, 0, 0) };
    disp[2] = parentHandle(earthNew[0], moonNew[0]);
    velo[2] = parentHandle(earthNew[1], moonNew[1]);

    std::cout << disp[2][0] << ',' << disp[2][1] << ',' << disp[2][2] << '\n' << velo[2][0] << ',' << velo[2][1] << ',' << velo[2][2] << '\n';

    std::ofstream myfile;
    myfile.open("output.csv"); // open a csv file in which to save the test simulation data

    myfile << "Time,Body,Disp (x),Disp (y),Disp (z),Velo (x),Velo (y),Velo (z),Acce (x),Acce (y),Acce (z),L Total,L,E,Bound,Ecce\n"; // set the header of the csv

    progressbar bar(3160); // creates progress bar for simulation 

    double totL{ 0 };
    double totE{ 0 };

    for (int its = 0; its < 3.16e+6; ++its) // number of time steps over which to run the simulation for
    {
        // calling simulation functions to perform the leapfrog integration
        updateAcce(disp, acce, mass, bodies);
        updateDisp(disp, velo, timestep, bodies);
        updateAcce(disp, acce, mass, bodies);
        updateVelo(velo, acce, timestep, bodies);
        updateDisp(disp, velo, timestep, bodies);

        time = time + timestep; // advance time each timestep

        if (fmod(time, 100000) == 0) // only output to the csv every 1000 iterations
        {
            for (int i = 0; i < bodies; ++i) // loop through each body
            {
                myfile << time << ',' << name[i]; // output the current time and body name to csv

                for (int j = 0; j < 3; ++j) // loop through each dimension
                {
                    myfile << ',' << disp[i][j] - disp[0][j]; // outputs the relative displacement of each body to the Sun
                }

                for (int j = 0; j < 3; ++j)
                {
                    myfile << ',' << velo[i][j] - velo[0][j]; // outputs the relative velocity of each body to the Sun
                }

                for (int j = 0; j < 3; ++j)
                {
                    myfile << ',' << acce[i][j] - acce[0][j]; // outputs the relative acceleration of each body to the Sun
                }

                totE = totalEner(disp, velo, mass, i, parent[i]);
                totL = totalMome(disp, velo, mass, i, parent[i]);

                myfile << ',' << angMome(disp, velo, mass, bodies) << ',' << totL << ',' << totE << ',' << bound(totE) << ',' << eccentricity(mass, totL, totE, i, parent[i]);

                myfile << '\n';
            }

            bar.update(); // update progress bar
        }

    }
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
