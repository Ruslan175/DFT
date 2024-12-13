#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

#define  EPS_VALUE      0.0001

using namespace std;

struct tDFTData
{
    double freq;
    double ampl;
    double re;
    double im;

    tDFTData(double f, double a, double r, double i) :  freq(f), ampl(a), re(r), im(i) {};
};


typedef     std::vector<double>     tSoundData;
typedef     std::vector<tDFTData>   tDftBuffer;


void getDFT(const tSoundData &data, int k, double &re, double &im);
double useDFT(const tSoundData &data, double dt, tDftBuffer &dft);
double func(double t);
void saveData(const char *name, tSoundData &data_t, tSoundData &data_f);
void useIDFT(const tDftBuffer &dft, tSoundData &data);


int main()
{
    tSoundData data_f;
    tSoundData data_t;
    tDftBuffer data_dft;
    const double dt = 0.01;
    const double T = 10.0 + 0.5 * dt;

    for (double t = 0.0; t < T; t += dt)
    {
        double f = func(t);
        if (abs(f) < EPS_VALUE)  f = 0.0;
        data_f.push_back(f);
        data_t.push_back(t);
    }
    const double max_amp = useDFT(data_f, dt, data_dft);

    // Show only main frequencies & updating
    for (auto &i : data_dft)
    {
        if (i.ampl/max_amp > 0.05)
        {
            cout << "f=" << i.freq << ":\tampl=" << i.ampl << ",\tre=" << i.re << ",\tim=" << i.im
                 << "\t(" <<  static_cast<int>(100.0 * i.ampl/max_amp) << "%)" << endl;
        }
        //if (i.ampl > 500.0)
        //{// Filter some frequencies
        //   i.re =  i.im = 0.0;
        //}
    }

    // Inverse DFT
    tSoundData data_idft;
    useIDFT(data_dft, data_idft);

    saveData("func.txt", data_t, data_f);
    saveData("func_inv.txt", data_t, data_idft);
    return 0;
}


double func(double t)
{
    return 2.0 * sin(3.0 * 2.0 * M_PI * t) - 0.80 * sin(2.0 * 2.0 * M_PI * t);
}


// DFT: https://www.sanfoundry.com/cpp-program-compute-discrete-fourier-transform-using-naive-approach/
// https://en.wikipedia.org/wiki/Discrete_Fourier_transform
void getDFT(const tSoundData &data, int k, double &re, double &im)
{
    const unsigned int sz = data.size();
    const double r = (2.0 * k * M_PI) / sz;
    re = 0.0;
    im = 0.0;

    for (unsigned int i = 0u; i < sz; ++i)
    {
        re += data[i] * cos(i * r);
        im += data[i] * sin(i * r);
    }
}


double useDFT(const tSoundData &data, double dt, tDftBuffer &dft)
{
    const double r = 1.0 / (dt * data.size());
    double max_amp = 0.0;
    double re, im;
    dft.clear();
    for (unsigned int k = 0u; k < data.size(); ++k)
    {
        const double freq = k * r;
        getDFT(data, k, re, im);
        const double amp = sqrt(re * re + im * im);
        
        if (amp > max_amp)
        {
            max_amp = amp;
        }

        dft.push_back(tDFTData(freq, amp, re, im));
    }
    return max_amp;
}


// File format is related to ...
// https://www.emathhelp.net/calculators/calculus-1/online-graphing-calculator/
void saveData(const char *name, tSoundData &data_t, tSoundData &data_f)
{
    ofstream file(name);
    if (false == file.is_open()) return;
    for (unsigned int i = 0; i < data_t.size(); ++i) file << data_t[i] << ", " << data_f[i] << "\n";
    file.close();
}


// https://www.geeksforgeeks.org/discrete-fourier-transform-and-its-inverse-using-c/
void useIDFT(const tDftBuffer &dft, tSoundData &data)
{
    data.clear();
    const unsigned int sz = dft.size();
    const double r = 2.0 * M_PI / sz;
    for (unsigned int n = 0u; n < sz; ++n)
    {   
        double val = 0.0;
        for (unsigned int k = 0u; k < sz; ++k) 
        {
            const double theta = r * k * n;
            val += (dft[k].re * cos(theta) + dft[k].im * sin(theta));
        }
        val /= sz;
        if (abs(val) < EPS_VALUE)  val = 0.0;
        data.push_back(val);
    }
}

