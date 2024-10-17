#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>

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

int main() {
    double precision = 1e-10;

    int monte_carlo_iterations = 0;
    int leibniz_iterations = 0;
    int bbpi_iterations = 0;

    std::cout << std::fixed << std::setprecision(10);
    std::cout << "Monte Carlo: " << calculate_pi_monte_carlo(precision, monte_carlo_iterations)
        << " (iterations: " << monte_carlo_iterations << ")" << std::endl;
    std::cout << "Leibniz: " << calculate_pi_leibniz(precision, leibniz_iterations)
        << " (iterations: " << leibniz_iterations << ")" << std::endl;
    std::cout << "Bailey-Borwein-Plouffe: " << calculate_pi_bbpi(precision, bbpi_iterations)
        << " (iterations: " << bbpi_iterations << ")" << std::endl;

    return 0;
}
