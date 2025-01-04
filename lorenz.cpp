//............................................. Start del código .............................................................................
/*
Este programa resuelve el sistema de Lorenz utilizando el método de Runge-Kutta de 4º orden.
Está diseñado para ser portátil y ejecutable en cualquier computadora con un compilador C++ estándar, es solo una implmentación simple.
El código está estructurado en clases para el sistema de Lorenz y el obvio solucionador de Runge-Kutta,
está la función principal para orquestar la simulación y generar los resultados en un archivo CSV.

This program solves the Lorenz system using a 4th order Runge-Kutta method.
It is designed to be portable and runnable on any computer with a standard C++ compiler.
The code is structured into classes for the Lorenz system and the Runge-Kutta solver,
with a main function to orchestrate the simulation and output the results to a CSV file.

By Felipe Correa Rodríguez.
*/
// ............................................. ............................................. .............................................

#include <vector>
#include <functional>
#include <iostream>
#include <fstream>
#include <stdexcept>

// .... Operaciones vectoriales ....
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> resultado(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        resultado[i] = a[i] + b[i];
    }
    return resultado;
}

std::vector<double> operator*(const std::vector<double>& a, double b) {
    std::vector<double> resultado(a.size());
    for (size_t i = 0; i < a.size(); ++i) {
        resultado[i] = a[i] * b;
    }
    return resultado;
}

std::vector<double> operator*(double b, const std::vector<double>& a) {
    return a * b;
}

std::vector<double>& operator+=(std::vector<double>& a, const std::vector<double>& b) {
    for (size_t i = 0; i < a.size(); ++i) {
        a[i] += b[i];
    }
    return a;
}

// .... Sistema de Lorenz .....
class sistema_lorenz {
public:
    sistema_lorenz(double sigma, double rho, double beta) : sigma_(sigma), rho_(rho), beta_(beta) {}

    void operator()(const std::vector<double>& estado, std::vector<double>& derivadas, double t) const {
        derivadas[0] = sigma_ * (estado[1] - estado[0]);
        derivadas[1] = estado[0] * (rho_ - estado[2]) - estado[1];
        derivadas[2] = estado[0] * estado[1] - beta_ * estado[2];
    }

private:
    double sigma_, rho_, beta_;
};

// .... Solucionador de Runge-Kutta .....
class solucionador_runge_kutta {
public:
    solucionador_runge_kutta(std::function<void(const std::vector<double>&, std::vector<double>&, double)> sistema_ode)
        : sistema_ode_(sistema_ode) {}

    void establecer_condiciones_iniciales(const std::vector<double>& estado_inicial) {
        estado_actual_ = estado_inicial;
    }

    void establecer_archivo_salida(const std::string& nombre_archivo) {
        archivo_salida_.open(nombre_archivo);
        if (!archivo_salida_.is_open()) {
            throw std::runtime_error("error al abrir el archivo de salida.");
        }
    }

    void integrar(double t_inicio, double t_fin, double dt) {
        double t = t_inicio;
        while (t < t_fin) {
            std::vector<double> k1(estado_actual_.size());
            sistema_ode_(estado_actual_, k1, t);
            std::vector<double> k2(estado_actual_.size());
            std::vector<double> estado_temporal = estado_actual_ + (dt / 2.0) * k1;
            sistema_ode_(estado_temporal, k2, t + dt / 2.0);
            std::vector<double> k3(estado_actual_.size());
            estado_temporal = estado_actual_ + (dt / 2.0) * k2;
            sistema_ode_(estado_temporal, k3, t + dt / 2.0);
            std::vector<double> k4(estado_actual_.size());
            estado_temporal = estado_actual_ + dt * k3;
            sistema_ode_(estado_temporal, k4, t + dt);
            estado_actual_ += (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
            t += dt;
            // salida del estado
            archivo_salida_ << t << ",";
            for (size_t i = 0; i < estado_actual_.size(); ++i) {
                archivo_salida_ << estado_actual_[i];
                if (i < estado_actual_.size() - 1) archivo_salida_ << ",";
            }
            archivo_salida_ << std::endl;
        }
    }

private:
    std::function<void(const std::vector<double>&, std::vector<double>&, double)> sistema_ode_;
    std::vector<double> estado_actual_;
    std::ofstream archivo_salida_;
};

// función principal
int main() {
    try {
        // definir condiciones iniciales
        std::vector<double> estado_inicial = {1.0, 1.0, 1.0};

        // definir parámetros para el sistema de lorenz
        double sigma = 10.0;
        double rho = 28.0;
        double beta = 8.0 / 3.0;

        // crear el sistema ode
        sistema_lorenz sistema_ode(sigma, rho, beta);

        // crear el solucionador
        solucionador_runge_kutta solucionador(sistema_ode);

        // establecer condiciones iniciales
        solucionador.establecer_condiciones_iniciales(estado_inicial);

        // establecer archivo de salida
        solucionador.establecer_archivo_salida("lorenz_output.csv");

        // integrar desde t_inicio hasta t_fin con paso de tiempo dt
        double t_inicio = 0.0;
        double t_fin = 50.0;
        double dt = 0.01;

        solucionador.integrar(t_inicio, t_fin, dt);
    } catch (const std::exception& e) {
        std::cerr << "se produjo una excepción: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}
// ............................................. ............................................. .............................................
// @FECORO
