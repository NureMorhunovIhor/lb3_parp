#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <windows.h>
#include <omp.h>

void task1() {
#ifdef _OPENMP
    printf("_OPENMP Defined\n");
#else
    printf("_OPENMP UnDefined\n");
#endif
}

void task2() {
    SYSTEM_INFO si;
    GetSystemInfo(&si);
    printf("Windows CPU cores: %u\n", si.dwNumberOfProcessors);

    int omp_threads = omp_get_max_threads();
    printf("OpenMP threads: %d\n", omp_threads);
}

double calculate_pi_monte_carlo(double precision, int& iterations) {
    int inside_circle = 0;
    int total_points = 0;
    double pi_estimate = 0.0;

    std::srand(std::time(0));

    while (true) {
        double x = static_cast<double>(std::rand()) / RAND_MAX * 2 - 1;
        double y = static_cast<double>(std::rand()) / RAND_MAX * 2 - 1;
        double distance = x * x + y * y;

        if (distance <= 1) {
            inside_circle++;
        }
        total_points++;
        pi_estimate = (static_cast<double>(inside_circle) / total_points) * 4;

        iterations++;

        if (std::abs(pi_estimate - 3.141592653589793) < precision) {
            return pi_estimate;
        }
    }
}

double calculate_pi_leibniz(double precision, int& iterations) {
    double pi_estimate = 0.0;
    int k = 0;

    while (true) {
        pi_estimate += (k % 2 == 0 ? 1.0 : -1.0) / (2 * k + 1);
        k++;
        iterations++;

        double current_pi = pi_estimate * 4;

        if (std::abs(current_pi - 3.141592653589793) < precision) {
            return current_pi;
        }
    }
}

double calculate_pi_bbpi(double precision, int& iterations) {
    double pi_estimate = 0.0;
    int k = 0;

    while (true) {
        double term = (1.0 / std::pow(16, k)) * (4.0 / (8 * k + 1) - 2.0 / (8 * k + 4) - 1.0 / (8 * k + 5) - 1.0 / (8 * k + 6));
        pi_estimate += term;
        iterations++;

        if (std::abs(pi_estimate - 3.141592653589793) < precision) {
            return pi_estimate;
        }
        k++;
    }
}

void task3_4_5() {
    double precision = 1e-10;

    int monte_carlo_iterations = 0;
    int leibniz_iterations = 0;
    int bbpi_iterations = 0;

    LARGE_INTEGER frequency, start, end;
    QueryPerformanceFrequency(&frequency);

    QueryPerformanceCounter(&start);
    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Monte Carlo: " << calculate_pi_monte_carlo(precision, monte_carlo_iterations)
        << " (iterations: " << monte_carlo_iterations << ")" << std::endl;
    QueryPerformanceCounter(&end);
    double monte_carlo_time = static_cast<double>(end.QuadPart - start.QuadPart) / frequency.QuadPart;

    QueryPerformanceCounter(&start);
    std::cout << "Leibniz: " << calculate_pi_leibniz(precision, leibniz_iterations)
        << " (iterations: " << leibniz_iterations << ")" << std::endl;
    QueryPerformanceCounter(&end);
    double leibniz_time = static_cast<double>(end.QuadPart - start.QuadPart) / frequency.QuadPart;

    QueryPerformanceCounter(&start);
    std::cout << "Bailey-Borwein-Plouffe: " << calculate_pi_bbpi(precision, bbpi_iterations)
        << " (iterations: " << bbpi_iterations << ")" << std::endl;
    QueryPerformanceCounter(&end);
    double bbpi_time = static_cast<double>(end.QuadPart - start.QuadPart) / frequency.QuadPart;

    std::cout << "\nExecution time for sequential methods:\n";
    std::cout << "Monte Carlo time: " << monte_carlo_time << " seconds\n";
    std::cout << "Leibniz time: " << leibniz_time << " seconds\n";
    std::cout << "Bailey-Borwein-Plouffe time: " << bbpi_time << " seconds\n";

    std::cout << "\nBest method based on iterations:\n";
    if (monte_carlo_iterations < leibniz_iterations && monte_carlo_iterations < bbpi_iterations) {
        std::cout << "Monte Carlo is the best method.\n";
    }
    else if (leibniz_iterations < monte_carlo_iterations && leibniz_iterations < bbpi_iterations) {
        std::cout << "Leibniz is the best method.\n";
    }
    else {
        std::cout << "Bailey-Borwein-Plouffe is the best method.\n";
    }
}

struct ThreadData {
    double precision;
    double pi_estimate;
    int iterations;
};

DWORD WINAPI calculate_pi_bbpi_parallel(LPVOID param) {
    ThreadData* data = (ThreadData*)param;
    double precision = data->precision;
    data->pi_estimate = 0.0;
    data->iterations = 0;
    int k = 0;

    while (true) {
        double term = (1.0 / std::pow(16, k)) * (4.0 / (8 * k + 1) - 2.0 / (8 * k + 4) - 1.0 / (8 * k + 5) - 1.0 / (8 * k + 6));
        data->pi_estimate += term;
        data->iterations++;

        if (std::abs(data->pi_estimate - 3.141592653589793) < precision) {
            return 0;
        }
        k++;
    }
}

void task6() {
    double precision = 1e-10;
    const int num_threads = 4;
    HANDLE threads[num_threads];
    ThreadData thread_data[num_threads];

    LARGE_INTEGER frequency, start, end;
    QueryPerformanceFrequency(&frequency);

    QueryPerformanceCounter(&start);

    for (int i = 0; i < num_threads; i++) {
        thread_data[i].precision = precision;
        threads[i] = CreateThread(NULL, 0, calculate_pi_bbpi_parallel, &thread_data[i], 0, NULL);
    }

    for (int i = 0; i < num_threads; i++) {
        WaitForSingleObject(threads[i], INFINITE);
        CloseHandle(threads[i]);
    }

    double pi_parallel = 0.0;
    int total_iterations = 0;

    for (int i = 0; i < num_threads; i++) {
        pi_parallel += thread_data[i].pi_estimate; 
        total_iterations += thread_data[i].iterations; 
    }

    pi_parallel /= num_threads; 

    QueryPerformanceCounter(&end);
    double parallel_time = static_cast<double>(end.QuadPart - start.QuadPart) / frequency.QuadPart;

    std::cout << "\nBailey-Borwein-Plouffe (parallel): " << std::fixed << std::setprecision(10) << pi_parallel << std::endl;
    std::cout << "Total iterations in parallel: " << total_iterations << std::endl;
    std::cout << "Parallel execution time: " << parallel_time << " seconds\n";
}

int main() {
    task1();
    task2();
    task3_4_5();
    task6();

    return 0;
}
