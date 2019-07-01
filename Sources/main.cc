#include <vector>
#include <iostream>

#include "ailibrary/ai.hh"
#include "planar3D.hh"

/*!
 \brief Версия программы
 \details Текущая версия программы (SemVer)
*/
std::string version("3.9.1");

#if defined(_MSC_VER)
    std::string compiler(ai::string("msc") + ai::string(_MSC_VER));
#elif defined(_COMPILER_)
    std::string compiler(TO_STRING(_COMPILER_));
#else
    std::string compiler("%compiler%");
#endif
#if defined(_TARGET_)
    std::string target(TO_STRING(_TARGET_));
#else
    std::string target("%target%");
#endif
#if defined(_TIMESTAMP_)
    std::string timestamp(TO_STRING(_TIMESTAMP_));
#else
    std::string timestamp("%timestamp%");
#endif
#if defined(_PLATFORM_)
    std::string buildVersion(
        version + ai::string("+") + TO_STRING(_PLATFORM_)
    );
#else
    std::string buildVersion(version + ai::string("+%platform%"));
#endif

#if defined(OPENMP)
    #include <omp.h>
#else
	//std::size_t omp_get_max_threads() { return 1; }
	int omp_get_max_threads() { return 1; }
#endif

/// \todo Defines: DEBUG, MEASURE_TIME, BUILD_DLL

/*!
 \brief Main-функция
 \details Первичная функция, считывающая параметры, переданные при запуске
 программы, и запускающая расчёт распространения трещины по планарной
 модели
*/
int main(const int argc, const char *argv[]){
    #if defined(BUILD_DLL)
    //std::cout.setstate(std::ios::failbit);
    DLL_Param DLL_Parametrs;

    std::vector< std::vector<double> > layers;
    std::vector< std::vector<double> > injection;
    std::string IdDesign;
    std::string IdStage;
    double modelingTime = -1;
    double Emit_time = -1;
    double Z_coordinate;
    std::string JSONstring;
    bool initialized = false;
    std::string strJson;



    ai::parseFileIntoString("init.json", strJson);

    DLL_Parametrs.J_String = strJson.c_str();

    std::cout << "start parce" << std::endl;
    ImportJSON(
        DLL_Parametrs,
        layers,
        injection
    );

    ai::saveMatrix("layers_in", layers);
    std::cout << "end parce" << std::endl;

    ApproximateLayers(layers, DLL_Parametrs.dx);  //перерасчет слоев и интерполирование на планаровскую сетку

    ai::saveMatrix("layers_appr", layers);
    // Пересчитываем значения для слоёв в величины СИ
    for (std::size_t i = 0; i < layers.size(); ++i) {
        layers[i][2] *= std::pow(10, 6);
        layers[i][3] *= std::pow(10, 9);
        layers[i][5] *= std::pow(10, -6);
    }


    // ai::saveMatrix("layers_planar", layers);
    //
    ai::saveMatrix("injection", injection);

    std::cout << "end init" << std::endl;

    std::cout << "start planar" << std::endl;
    int a = planar3D(
        DLL_Parametrs,
        layers,
        injection,
        initialized
    );

    return 0;
    #endif

    bool override = false;
    bool saveSteps = false;
    bool considerElasticModulusContrast = false;

    std::string pathToLayersFile("layers.txt");
    std::string pathToInjectionFile("injection.txt");
    std::string pathToImportFolder = std::string("./InitialConditions");

    std::size_t numberOfThreads = (std::size_t) omp_get_max_threads();

    int meshScalingCounter = -1;

    double cellSize = -1.;
    double meshSize = 30.;
    double time = 25.0;
    double timeStep = -1.;
    double timeScale = 60.;

    epsilon = std::pow(10., -8);
    regime = VISCOSITY;

    for(int i = 1; i < argc; ++i){
        if(
            std::string("-v") == std::string(argv[i])
            || std::string("--version") == std::string(argv[i])
        ){
            std::cout << "  Build: "  << buildVersion << "." << std::endl;
            std::cout << "  Compiler: "  << compiler << "." << std::endl;
            std::cout << "  Target: "  << target << "." << std::endl;
            std::cout << "  AiLibrary: " << ai::getVersion() << "."
                << std::endl;
            std::cout << "  Compilation timestamp: "  << timestamp << "."
                << std::endl;
            std::cout << "  Additional features:" << std::endl;
            std::cout << "\tIMEX - ";
            #if defined(IMEX)
            std::cout << "yes" << std::endl;
            #else
            std::cout << "no" << std::endl;
            #endif
            std::cout << "\tDEBUG - ";
            #if defined(DEBUG)
            std::cout << "yes" << std::endl;
            #else
            std::cout << "no" << std::endl;
            #endif
            std::cout << "\tEIGEN - ";
            #if defined(USE_EIGEN)
            std::cout << "yes" << std::endl;
            #else
            std::cout << "no" << std::endl;
            #endif
            std::cout << "\tFLUID MARKERS - ";
            #if defined(FLUID_MARKERS)
            std::cout << "yes" << std::endl;
            #else
            std::cout << "no" << std::endl;
            #endif
            std::cout << "\tPROPPANT MARKERS - ";
            #if defined(PROPPANT_MARKERS)
            std::cout << "yes" << std::endl;
            #else
            std::cout << "no" << std::endl;
            #endif
            std::cout << "\tTRUE ELASTIC CONTRAST - ";
            #if defined(TRUE_ELASTIC_CONTRAST)
            std::cout << "yes" << std::endl;
            #else
            std::cout << "no" << std::endl;
            #endif

            return 0;
        }

        if(std::string("--name") == std::string(argv[i])){
            std::cout << "Planar3D-LL" << std::endl;

            return 0;
        }

        if(
            std::string("-h") == std::string(argv[i])
            || std::string("--help") == std::string(argv[i])
        ){
            std::cout << "usage: planar3D [options]"
                << std::endl
                << "    -h  --help            print this usage and exit"
                << std::endl
                << "    -v  --version         print build info and exit"
                << std::endl
                << "    --name                print program name and exit"
                << std::endl
                << "    --list-errors         print possible errors ans exit"
                << std::endl << std::endl

                << "  Time parameters" << std::endl
                << "    --time=<value>        modeling time [double, min]"
                << std::endl
                << "    --time-step=<value>   time step [double, min]"
                << std::endl
                << "    --time-scale=<value>  time scale [double, s]"
                << std::endl << std::endl

                << "  Mesh parameters" << std::endl
                << "    --mesh-size=<value>   mesh size in cells [uint, n/d]"
                << std::endl
                << "    --cell-size=<value>   cell size in meters [double, m]"
                << std::endl
                << "    --mesh-scale=<value>  how many times mesh scaling is "
                << "allowed [uint, n/d]"
                << std::endl
                << "    --mesh-scale          allow infinite mesh scaling"
                << std::endl
                << std::endl << std::endl

                << "  Initial fracture regime" << std::endl
                << "    --viscosity           viscosity dominated regime "
                << "{default}"
                << std::endl
                << "    --toughness           toughness dominated regime "
                << std::endl
                << "    --leak-off            leak-off dominated regime "
                << std::endl << std::endl

                << "  Initial conditions" << std::endl
                << "    --layers=<path>       path to a txt-file with a layer "
                << "description [string]"
                << std::endl
                << "    --injection=<path>    path to a txt-file with a plan "
                << "of injections [string]"
                << std::endl << std::endl

                << "  Other flags" << std::endl
                << "    --threads=<value>     number of parallel openmp"
                << " threads [uint, n/d]"
                << std::endl
                << "    --elastic-contrast    consider the effect of elastic "
                << "modulus contrast"
                << std::endl
                << "    --save-steps          save concentration and opening "
                << "at every major time step"
                << std::endl
                << "    --override            ignore possible warnings"
                << std::endl
                << "    --import=<path>       import data from folder [string]"
                << std::endl
                << "    --import              default import from ./Results"
                << std::endl;

            return 0;
        }

        if("--list-errors" == std::string(argv[i])){
            std::cout << "  User input errors" << std::endl
                << "    Code 11. Cell size is less than 1 meter."
                << std::endl
                << "    Code 12. Mesh size is less than 10 cells."
                << std::endl
                << "    Code 13. Mesh size is not a positive integer."
                << std::endl
                << "    Code 17. Time step isn't a positive value."
                << std::endl << std::endl

                << "  Enviroment errors" << std::endl
                << "    Code 21. Cannot find %folderName%."
                << std::endl
                << "    Code 22. Cannot open %fileName%."
                << std::endl
                << "    Code 23. Format error in %fileName%."
                << std::endl << std::endl

                << "  Calculating errors" << std::endl
                << "    Code 31. Calculations failed due to unmet "
                << "mathematical condition."
                << std::endl
                << "    Code 32. Automodel solution cannot be applied. "
                << std::endl;

            return 0;
        }

        if(
            ai::assignAbsDoubleParameter(argv[i], "--time=", time)
            || ai::assignAbsDoubleParameter(argv[i], "--time-step=", timeStep)
            || ai::assignAbsDoubleParameter(argv[i], "--time-scale=", timeScale)
            || ai::assignAbsDoubleParameter(argv[i], "--cell-size=", cellSize)
            || ai::assignAbsDoubleParameter(argv[i], "--mesh-size=", meshSize)
            || ai::assignStringParameter(
                argv[i],
                "--layers=",
                pathToLayersFile
            )
            || ai::assignStringParameter(
                argv[i],
                "--injection=",
                pathToInjectionFile
            )
            || ai::assignParameter(
                argv[i],
                "--mesh-scale=",
                meshScalingCounter
            )
            || ai::assignParameter(
                argv[i],
                "--import=",
                pathToImportFolder
            )
            || ai::assignParameter(argv[i], "--threads=", numberOfThreads)
        ){
            continue;
        }

        if("--viscosity" == std::string(argv[i])){
            regime = VISCOSITY;

            continue;
        }

        if("--toughness" == std::string(argv[i])){
            regime = TOUGHNESS;

            continue;
        }

        if("--leak-off" == std::string(argv[i])){
            regime = LEAK_OFF;

            continue;
        }

        if("--mesh-scale" == std::string(argv[i])){
            meshScalingCounter = 1000;

            continue;
        }

        if("--import" == std::string(argv[i])){
            pathToImportFolder = "./Results";

            continue;
        }

        if("--elastic-contrast" == std::string(argv[i])){
            considerElasticModulusContrast = true;

            continue;
        }

        if("--save-steps" == std::string(argv[i])){
            saveSteps = true;

            continue;
        }

        if("--override" == std::string(argv[i])){
            override = true;

            continue;
        }
    }


    if(1. > cellSize && -1. != cellSize && !override){
        std::cerr << "Cell size should be at least 1 meter" << std::endl;

        return 11;
    }

    if(10. > meshSize && !override){
        std::cerr << "Mesh size should be at least 10 cells." << std::endl;

        return 12;
    }

    if(meshSize != (int) meshSize && !override){
        std::cerr << "Mesh size should be a positive integer." << std::endl;

        return 13;
    }

    if(0 > timeStep && -1. != timeStep && !override){
        std::cerr << "Time step should be a positive value." << std::endl;

        return 13;
    }

    if(0 != (int) meshSize % 2){
        meshSize = round(meshSize + 1);
    }

    /// \todo time, mesh, data, settings
    #if defined(BUILD_DLL)
    return 0;
    #else
    return planar3D(
       time, timeStep, timeScale, cellSize, meshSize, meshScalingCounter,
       pathToLayersFile, pathToInjectionFile, pathToImportFolder,
       considerElasticModulusContrast, saveSteps, false,
       numberOfThreads, override
    );
    #endif
}
