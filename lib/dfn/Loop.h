#ifndef LOOP_H
#define LOOP_H

#include <mechsys/dfn/Domain.h>
#include <omp.h>

namespace DFN
{
    class Loop
    {
    public:
        double times;                   ///< loop times, each time DFN modeling, the density will be increased compared to last time DFN modeling
        double RatioLR;                 ///< ratio of L and Ra, L is model edge size, Ra is the average radius of min and max ourter circle radius
        double R_a;                     ///< the average radius of min and max ourter circle radius
        double R_low;                   ///< minimum radius of outer circle radius
        double R_up;                    ///< maximum radius of outer circle radius
        double L;                       ///< model size
        size_t Nproc;                   ///< num of threads
        size_t nt;                      ///< should not larger than times! which model is the model shown to us
        size_t nk;                      ///< when nk DFN models are finished, output one time
        size_t nv_MC_TIMES;             ///< each density, the MC times
        double nx;                      ///< the increment of fracture number regard to each DFN modeling time
        double Percolation_parameter_c; ///< when percolation probability equals to 0.5, the percolation parameter is
        Array<double> P32_total_1;      ///< the next several 1D arrays store outputs
        Array<double> P32_connected_1;
        Array<double> P30_1;
        Array<double> Percolation_parameter_1;
        Array<double> Ratio_of_P32_1;
        Array<double> Percolation_probability_1;

    public:
        void Loop_create_DFNs(double vari, double random_seed, String str_ori, String percolation_direction, double array13[7]);
        void Data_output_stepBYstep(size_t times, char const *FileKey, double alpha, double L, double increment_fracture, double R_low, double R_up, double P32_total_B, double P32_connected_B, double P30_B, double Percolation_parameter_B, double Ratio_of_P32_B, double Percolation_probability_B);
        void Sign_of_finding_pc(char const *FileKey); ///< if Pc is found, outputs a file
    };

    void Loop::Loop_create_DFNs(double vari, double random_seed, String str_ori, String percolation_direction, double array13[7])
    {

        //times = 2000;                     // times times DFN, each DFN time has nv times MC simulation, each DFN time has the same density, which is larger than the last DFN time
        //size_t nt = 200;                  //should not larger than times! which model is the model shown to us
        //size_t nk = 1;                    //when nk DFN models are finished, output one time
        size_t nv = nv_MC_TIMES; //each density, the MC times
        //double nx = 10;                   //the increment of fracture number regard to each DFN modeling time
        double alpha_for_powerLaw = vari; // exponent for power law

        double alpha_radius = -alpha_for_powerLaw; //this is used to generate random numbers
        //RatioLR = 10;

        //R_a = 6; //it is not the expectation, actually it is the average of lower and upper boundaries of fracture radius

        //R_low = R_a / 2;
        //R_up = R_a * 1.5;
        L = RatioLR * R_a;

        double array11[3][2] = {-L * 0.5 - R_up, L * 0.5 + R_up, -L * 0.5 - R_up, L * 0.5 + R_up, -L * 0.5 - R_up, L * 0.5 + R_up}; //min x, max x; min y, max y; min z, max z; (they are DOMAIN size, not MODEL size!!!)
        double array12[4] = {0, alpha_radius, R_low, R_up};
        //   double array13[7] = {30, 20, 20, 10, 100.1, 0.5, 30.1};																				//mean radius(diagonal line), standard error, min radius, max radius
        double model_size[6] = {-L * 0.5, L * 0.5, -L * 0.5, L * 0.5, -L * 0.5, L * 0.5}; // MODEL size, not DOMAIN
        //double model_side_length = model_size[1] - model_size[0];
        //double outer_circle_radius = array12[0];

        size_t np = 0;
        size_t njk = 0;

        while (np < times)
        {
            np++;
            size_t n = np * nx;

            double nf = 0;

            Array<double> P32_total_A;
            Array<double> P32_connected_A;
            Array<double> P30_A;
            Array<double> Percolation_parameter_A;
            Array<double> Ratio_of_P32_A;
            Array<double> Percolation_probability_A;

            P32_total_A.resize(nv);
            P32_connected_A.resize(nv);
            P30_A.resize(nv);
            Percolation_parameter_A.resize(nv);
            Ratio_of_P32_A.resize(nv);
            Percolation_probability_A.resize(nv);

#pragma omp parallel for schedule(static) num_threads(Nproc)
            for (size_t i = 0; i < nv; i++)
            {
                DFN::Domain dom;
                dom.Create_whole_model(n, random_seed, model_size, str_ori, array11, array12, array13); ///uniform means oritation data are generated uniformly, so, actually, array13 is input but not used

                size_t z = dom.Identify_percolation_clusters(percolation_direction);
                dom.Connectivity_uniform_orientation(percolation_direction);

                P32_total_A[i] = (dom.P32_total);
                P32_connected_A[i] = (dom.P32_connected);
                P30_A[i] = (dom.P30);
                Percolation_parameter_A[i] = (dom.Percolation_parameter);
                Ratio_of_P32_A[i] = (dom.Ratio_of_P32);

                if (z == 1)
                    Percolation_probability_A[i] = 1;
                else
                    Percolation_probability_A[i] = 0;

                if (np == nt && i == nv - 1)
                {
                    dom.WriteFrac("tdfn01");
                    dom.PlotMatlab_DFN("tdfn01_DFN");
                    dom.PlotMatlab_DFN_and_Intersection("tdfn01_DFN_and_Intersections");
                    dom.PlotMatlab_ORI_SCATTER("tdfn01_ORI_SCATTER");
                    //dom.PlotMatlab_Traces_on_Model_surfaces("tdfn01_Trace_on_surfaces");
                    dom.PlotMatlab_DFN_Highlight_Cluster("tdfn01_DFN_Highlight_Cluster");
                    dom.PLotMatlab_DFN_Cluster_along_a_direction("tdfn01_DFN_Z_clusters", "z");
                    dom.PlotMatlab_Radius_and_Area_kstest("tdfn01_DFN_Fracture_Radius_and_Area");
                    dom.PlotMatlab_Radius_and_Perimeter_kstest("tdfn01_DFN_Fracture_Radius_and_Perimeter");
                    dom.DataFile_Radius_AreaAndPerimeter("tdfn01_DFN_Radius_AreaAndPerimeter");
                }

                if (i == nv / 2 && np % nk == 0)
                {
                    using namespace std;
                    std::cout << "The Model NO." << np << " has been created! "
                              << "Times: " << i << "; Alpha: " << vari << "; thread: " << omp_get_thread_num() << std::endl;
                }
                if (i == nv - 1 && np % nk == 0)
                {
                    using namespace std;
                    std::cout << "The Model NO." << np << " has been created! "
                              << "Times: " << i << "; Alpha: " << vari << "; thread: " << omp_get_thread_num() << std::endl;
                }
            }

            double P32_total_B = 0;
            for (size_t i = 0; i < P32_total_A.size(); ++i)
            {
                P32_total_B = P32_total_B + P32_total_A[i];
            }
            P32_total_1.Push(P32_total_B / P32_total_A.size());

            double P32_connected_B = 0;
            for (size_t i = 0; i < P32_connected_A.size(); ++i)
            {
                P32_connected_B = P32_connected_B + P32_connected_A[i];
            }
            P32_connected_1.Push(P32_connected_B / P32_connected_A.size());

            double P30_B = 0;
            for (size_t i = 0; i < P30_A.size(); ++i)
            {
                P30_B = P30_B + P30_A[i];
            }
            P30_1.Push(P30_B / P30_A.size());

            double Percolation_parameter_B = 0;
            for (size_t i = 0; i < Percolation_parameter_A.size(); ++i)
            {
                Percolation_parameter_B = Percolation_parameter_B + Percolation_parameter_A[i];
            }
            Percolation_parameter_1.Push(Percolation_parameter_B / Percolation_parameter_A.size());

            double Ratio_of_P32_B = 0;
            for (size_t i = 0; i < Ratio_of_P32_A.size(); ++i)
            {
                Ratio_of_P32_B = Ratio_of_P32_B + Ratio_of_P32_A[i];
            }
            Ratio_of_P32_1.Push(Ratio_of_P32_B / Ratio_of_P32_A.size());

            for (size_t i = 0; i < Percolation_probability_A.size(); ++i)
            {
                if (Percolation_probability_A[i] == 1)
                    nf++;
            }
            Percolation_probability_1.Push(nf / ((double)(nv * 1.00)));
            double Percolation_probability_B = nf / ((double)(nv * 1.00));
            Data_output_stepBYstep(np, "tdfn01_datafile_step_by_step", alpha_for_powerLaw, L, nx, R_low, R_up, P32_total_B, P32_connected_B, P30_B, Percolation_parameter_B, Ratio_of_P32_B, Percolation_probability_B);

            if (njk == 0 && Percolation_probability_B > 0.49999)
            {
                njk++;
                Percolation_parameter_c = Percolation_parameter_B;
                std::cout << "\n**********Found Pc**********\n\n";
                Sign_of_finding_pc("Pc_Found");
            }

            if (njk != 0 && Percolation_parameter_B > 2 * Percolation_parameter_c)
            {
                std::cout << "\n**********Found two times Pc**********\n\n";
                break;
            }
        }
        PLotMatlab_DFN_Connectivity("tdfn01_Connectivity", Percolation_parameter_1, Percolation_probability_1, vari);

        std::cout << "Loop finished!\n";
        //Datafile_output("tdfn01_Datafile", model_side_length, P32_total_1, P32_connected_1, P30_1, Percolation_parameter_1, Ratio_of_P32_1, Percolation_probability_1);
    }; // namespace DFN

    inline void Loop::Data_output_stepBYstep(size_t times, char const *FileKey, double alpha, double L, double increment_fracture, double R_low, double R_up, double P32_total_B, double P32_connected_B, double P30_B, double Percolation_parameter_B, double Ratio_of_P32_B, double Percolation_probability_B)
    {
        if (times == 1)
        {
            std::ostringstream oss;
            oss << "Exponent:"
                << "\t" << alpha << "\n";
            oss << "Model_edge"
                << "\t" << L << "\n";
            oss << "Fracture_increment:"
                << "\t" << increment_fracture << "\n";
            oss << "Min_radius_of_outer_circle"
                << "\t" << R_low << "\n";
            oss << "Max_radius_of_outer_circle"
                << "\t" << R_up << "\n"
                << "\n";
            oss << "P32_total"
                << "\t"
                << "P32_connected"
                << "\t"
                << "P30"
                << "\t"
                << "Percolation_parameter"
                << "\t"
                << "Ratio_of_P32"
                << "\t"
                << "Percolation_probability"
                << "\n";
            oss << P32_total_B << "\t" << P32_connected_B << "\t" << P30_B << "\t" << Percolation_parameter_B << "\t" << Ratio_of_P32_B << "\t" << Percolation_probability_B << "\n";
            String fn(FileKey);
            fn.append(".txt");
            std::ofstream of(fn.CStr(), std::ios::out);
            of << oss.str();
            of.close();
        }
        else
        {
            std::ostringstream oss;
            oss << P32_total_B << "\t" << P32_connected_B << "\t" << P30_B << "\t" << Percolation_parameter_B << "\t" << Ratio_of_P32_B << "\t" << Percolation_probability_B << "\n";
            String fn(FileKey);
            fn.append(".txt");
            std::ofstream of(fn.CStr(), std::ios::app);
            of << oss.str();
            of.close();
        }
    };

    void Loop::Sign_of_finding_pc(char const *FileKey)
    {
        std::ostringstream oss;
        oss << "Pc has been found: ";
        oss << Percolation_parameter_c;

        String fn(FileKey);
        fn.append(".txt");
        std::ofstream of(fn.CStr(), std::ios::out);
        of << oss.str();
        of.close();
    }
}; // namespace DFN
#endif
