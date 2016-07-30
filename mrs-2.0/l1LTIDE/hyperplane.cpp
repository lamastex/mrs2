#include <iostream>
#include <sstream>

#include "DensityTree/CGALTypeDefs.hpp"
#include "DensityTree/PointUtils.hpp"
#include "DensityTree/bsp.hpp"


void test_bsp_tree(int dim, int num_points)
{
    //TODO: Make timings for various numbers of points and dimensions
    //      Use std::chrono
    //      std::cout << dimension << num points << time

    std::cout << "dim = " << dim << std::endl;
    std::cout << "num_points = " << num_points << std::endl;

    std::vector<Point> points = randPoints(dim, num_points, GeneratorType::Uniform);

    std::vector<Point*> point_pntr;
    for (auto&& p: points) {
        point_pntr.push_back(&p);
//        std::cout << p << std::endl;
    }

    BSP_Tree tree = BSP_Tree(point_pntr, 20);

    std::cout << std::endl;
    tree.root.print_tree(4, 1);
    int max_count = tree.root.Max_Num();
    int min_count = tree.root.Min_Num();
    double ratio = static_cast<double> (max_count) / static_cast<double> (min_count);
    std::cout << max_count << std::endl;
    std::cout << min_count << std::endl;
    std::cout << ratio << std::endl;
}

int main(int argc, char* argv[])
{
//    std::cout << "argc= " << argc << std::endl;

    int num_points;
    int dim;
    if (argc == 3) {

        std::istringstream arg1(argv[1]);
        std::istringstream arg2(argv[2]);
        if (!(arg1 >> num_points) || !(arg2 >> dim)) {
            throw "Invalid command line arguments";
        }
    } else {
        num_points = 10;
        dim = 2;
    }

    test_bsp_tree(dim, num_points);

    return 0;
}

