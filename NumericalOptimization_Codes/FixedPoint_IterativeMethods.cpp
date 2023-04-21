#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <cmath>
#include <tuple>
using namespace std;

vector<vector<double>> readfile_matrix(string filename)
{
    int m=2, n=2;
    // cout << "How many rows in array?";
    // cin >> m;
    // cout << "How many columns in array?";
    // cin >> n;
    vector<vector<double>> data;

    ifstream file(filename);
    if(file.is_open())
    {
        for(int i = 0; i < m; i++)
        {
            vector<double> vt2;
            double val;
            for(int j = 0; j < n; j++)
            {
                file >> val;
                vt2.push_back(val);
            }
            data.push_back(vt2);
        }
    }
    return data;
};

vector<double> readfile_vector_new(string filename)
{
    vector<double> data;
    int n = 2;
    double val;
    ifstream file(filename);
    if(file.is_open())
    {
        for(int j = 0; j < n; j++)
        {
            file >> val;
            data.push_back(val);
        }
    }
    return data;
};

vector<double> readFile_vector(string filename){
    vector<double> data;
    ifstream inFile(filename); // Data.txt
    if (inFile.is_open()){
        string line;
        while( getline(inFile,line))
        {
            stringstream ss(line);

            string num;
            while( getline(ss,num,',') ){
                data.push_back(stod(num.c_str()));
            };
        };
    };
return data;
};

class IterativeMethods
{
    private:
    // M^-1 A x = M^-1 b
    vector<vector<double>> A_matrix_dp;
    vector<vector<float>> A_matrix_sp;
    vector<vector<double>> Preconditioner_dp;
    vector<vector<float>> Preconditioner_sp;
    vector<double> b_vector_dp;
    vector<float> b_vector_sp;
    vector<double> solution_vector_dp;
    vector<float> solution_vector_sp_RF, solution_vector_sp_SD, solution_vector_sp_CG, solution_vector_sp_GS, solution_vector_sp_SD_slow, solution_vector_sp_CG_newton, solution_vector_sp_quasi_newt_BFGS;
    bool includePrecondition = false;
    vector<double> metric_known, metric_unknown, metric_precondition; 

    void solve_exact_via_steepest_descent()
    {
        solution_vector_dp = Iterate_dp("SD", 
        A_matrix_dp, 
        b_vector_dp, 
        solution_vector_dp,
        Preconditioner_dp,
        includePrecondition
        );
    };

    void solve_exact_via_conjugate_gradient()
    {
        solution_vector_dp = Iterate_dp("CG", 
        A_matrix_dp, 
        b_vector_dp, 
        solution_vector_dp,
        Preconditioner_dp,
        includePrecondition
        );
    };

    void solve_exact_via_richardson_gauss_seidel()
    {
        solution_vector_dp = Iterate_dp("RF", 
        A_matrix_dp, 
        b_vector_dp, 
        solution_vector_dp,
        Preconditioner_dp,
        includePrecondition
        );
    };

    void solve_via_steepest_descent()
    {
        solution_vector_sp_SD = Iterate("SD", 
        A_matrix_sp, 
        b_vector_sp, 
        solution_vector_sp_SD,
        Preconditioner_sp,
        includePrecondition
        );
    };

    void solve_via_gauss_southwell()
    {
        solution_vector_sp_GS = Iterate("GS", 
        A_matrix_sp, 
        b_vector_sp, 
        solution_vector_sp_GS,
        Preconditioner_sp,
        includePrecondition
        );
    };

    void solve_via_slow_steepest_descent()
    {
        solution_vector_sp_SD_slow = Iterate("SD_slow", 
        A_matrix_sp, 
        b_vector_sp, 
        solution_vector_sp_SD_slow,
        Preconditioner_sp,
        includePrecondition
        );
    };

    void solve_via_conjugate_gradient()
    {
        solution_vector_sp_CG = Iterate("CG", 
        A_matrix_sp, 
        b_vector_sp, 
        solution_vector_sp_CG,
        Preconditioner_sp,
        includePrecondition
        );
    };

    void solve_via_richardson_gauss_seidel()
    {
        solution_vector_sp_RF = Iterate("RF", 
        A_matrix_sp, 
        b_vector_sp, 
        solution_vector_sp_RF,
        Preconditioner_sp,
        includePrecondition
        );
    };

    void solve_via_newton_trunc_CG()
    {
        solution_vector_sp_CG_newton = Iterate("CG_newton", 
        A_matrix_sp, 
        b_vector_sp, 
        solution_vector_sp_CG_newton,
        Preconditioner_sp,
        includePrecondition
        );
    };

    void solve_via_quasi_newt_BFGS()
    {
        solution_vector_sp_quasi_newt_BFGS = Iterate("BFGS", 
        A_matrix_sp, 
        b_vector_sp, 
        solution_vector_sp_quasi_newt_BFGS,
        Preconditioner_sp,
        includePrecondition
        );
    };

    public:    

        // Overloaded Constructors
        // If Supplied Preconditioner
        IterativeMethods(vector<vector<double>> inputMatrix, vector<double> inputVector_b, vector<double> startVector, vector<vector<double>> preConditioner, bool includeConditioner)
        {
            includeConditioner = includeConditioner;
            copy(inputMatrix.begin(), inputMatrix.end(),back_inserter(A_matrix_dp));
            copy(preConditioner.begin(), preConditioner.end(),back_inserter(Preconditioner_dp));
            for (auto&& v : A_matrix_dp) A_matrix_sp.emplace_back(begin(v), end(v));
            for (auto&& v : Preconditioner_dp) Preconditioner_sp.emplace_back(begin(v), end(v));    
            b_vector_dp = inputVector_b;
            transform(b_vector_dp.begin(), b_vector_dp.end(),back_insert_iterator(b_vector_sp),[](double val) -> float{return (float)val;});
            solution_vector_dp = startVector;
            transform(solution_vector_dp.begin(), solution_vector_dp.end(),back_insert_iterator(solution_vector_sp_RF),[](double val) -> float{return (float)val;});
            solution_vector_sp_SD = solution_vector_sp_RF;
            solution_vector_sp_CG = solution_vector_sp_RF;
            solution_vector_sp_GS = solution_vector_sp_RF;
            solution_vector_sp_SD_slow = solution_vector_sp_RF;
            solution_vector_sp_CG_newton = solution_vector_sp_RF;
            solution_vector_sp_quasi_newt_BFGS = solution_vector_sp_RF;
        };
        // If NOT Supplied Preconditioner
        IterativeMethods(vector<vector<double>> inputMatrix, vector<double> inputVector_b, vector<double> startVector)
        {
            copy(inputMatrix.begin(), inputMatrix.end(),back_inserter(A_matrix_dp));
            A_matrix_sp.reserve(A_matrix_dp.size());
            for (auto&& v : A_matrix_dp) A_matrix_sp.emplace_back(begin(v), end(v));
            b_vector_dp = inputVector_b;
            transform(b_vector_dp.begin(), b_vector_dp.end(),back_insert_iterator(b_vector_sp),[](double val) -> float{return (float)val;});
            solution_vector_dp = startVector;
            transform(solution_vector_dp.begin(), solution_vector_dp.end(),back_insert_iterator(solution_vector_sp_RF),[](double val) -> float{return (float)val;});
            solution_vector_sp_SD = solution_vector_sp_RF;
            solution_vector_sp_CG = solution_vector_sp_RF;
            solution_vector_sp_GS = solution_vector_sp_RF;
            solution_vector_sp_SD_slow = solution_vector_sp_RF;
            solution_vector_sp_CG_newton = solution_vector_sp_RF;
            solution_vector_sp_quasi_newt_BFGS = solution_vector_sp_RF;
        };

        vector<double> Iterate_dp(
            string method, 
            vector<vector<double>> A_matrix, 
            vector<double> b_vector, 
            vector<double> x_vector, 
            vector<vector<double>> preConditioner,
            bool condition
            )
        {
            vector<double> vector_sol = x_vector;
            vector<vector<double>> Conditioner_factored = cholesky(preConditioner);
            vector<double> dummy_vector = matrix_vector_multiply(A_matrix, x_vector);
            vector<double> residual = b_vector;
            vector<double> r_prior = residual;
            transform(residual.begin(), residual.end(), dummy_vector.begin(), residual.begin(), minus<double>());
            int k = 0;
            vector<double> p = residual;

            if (method != "CG_newton")
            {
                bool zeros = all_of(residual.begin(), residual.end(), [](double i) { return i==0.0; });
                while (zeros == false && k<=50)
                {
                    p = get_direction(method, condition, k, Conditioner_factored, residual, r_prior, p);
                    vector<double> q = matrix_vector_multiply(A_matrix, p);
                    double alpha = inner_product(p, residual) / inner_product(p, q);
                    transform(p.begin(), p.end(), p.begin(), [alpha](double &c){ return c*alpha; }); // p = alpha * p
                    transform(vector_sol.begin(), vector_sol.end(), p.begin(), vector_sol.begin(), plus<double>()); // x = x + alpha*p
                    transform(q.begin(), q.end(), q.begin(), [alpha](double &c){ return c*alpha; }); // q = alpha * q
                    transform(residual.begin(), residual.end(), q.begin(), residual.begin(), minus<double>()); // r = r - alpha*q
                    k++;
                    zeros = all_of(residual.begin(), residual.end(), [](double i) { return i==0.0; });
                };

            }else if (method == "CG_newton")
            {
                bool zeros = all_of(residual.begin(), residual.end(), [](double i) { return i==0.0; });
                while  (zeros == false && k<=50)
                {
                
                    // Truncated Conjugate Gradient Inner Loop
                    bool terminate = false;
                    vector<double> d = residual;
                    int size = vector_sol.size();
                    vector<double> z(size);
                    fill(z.begin(), z.end(), 0);
                    int i = 0; 
                    while (!terminate) 
                    {
                        double sigma = inner_product(residual, residual);
                        vector<double> q = matrix_vector_multiply(A_matrix, d);
                        double mu = inner_product(d, q);
                        if (mu < 0)
                        {
                            p =  z;
                            break;
                        };
                        double alpha = sigma / mu;
                        transform(d.begin(), d.end(), d.begin(), [alpha](double &c){ return c*alpha; }); // d = alpha * d
                        transform(z.begin(), z.end(), d.begin(), z.begin(), plus<double>()); // z = z + alpha*d
                        transform(q.begin(), q.end(), q.begin(), [alpha](double &c){ return c*alpha; }); // q = alpha * q
                        transform(residual.begin(), residual.end(), q.begin(), residual.begin(), minus<double>()); // r = r - alpha*q
                        d = get_direction("CG", condition, i, Conditioner_factored, residual, r_prior, d);

                        double norm_prior = norm_two_dp(r_prior);
                        terminate = (norm_two_dp(residual) < min(0.5, sqrt(norm_prior)) * norm_prior);
                        i++;
                    };

                    float sig = inner_product(residual, residual);
                    vector<double> que = matrix_vector_multiply(A_matrix, p);
                    float mew = inner_product(p, que);
                    double alph = sig / mew;
                    transform(p.begin(), p.end(), p.begin(), [alph](double &c){ return c*alph; }); // d = alpha * d
                    transform(vector_sol.begin(), vector_sol.end(), p.begin(), vector_sol.begin(), plus<double>()); // x = x + alpha*p
                    k++;
                    zeros = all_of(residual.begin(), residual.end(), [](double i) { return i==0.0; });
                };

            };

            return vector_sol;
        };

        vector<vector<double>> cholesky(vector<vector<double>> matrix)
        {
            // Factor Matrix In Place
            int m_rows = (int)matrix.size();
            int n_columns = (int)matrix[0].size();

            for (int col = 0; col < n_columns; col++)
            {
                if (matrix[col][col] == 0)
                {
                    throw runtime_error("Math error: Attempted to divide by Zero. Consider Pivoting Matrix\n");
                }else
                {
                    // Root of Top Left
                    matrix[col][col] = sqrt(matrix[col][col]);
                    // Calculate Lower Trapeziodal Elements
                    for (int row = col+1; row < m_rows; row++)
                    {
                        matrix[row][col] = matrix[row][col] / matrix[col][col];
                    }
                    // Update SubMatrix Elements Before Next Iteration
                    for (int sub_row = col+1; sub_row < m_rows; sub_row++)
                    {
                        for (int sub_col = col+1; sub_col < n_columns; sub_col++)
                        {
                            matrix[sub_row][sub_col] = matrix[sub_row][sub_col] - ( matrix[sub_row][col] * matrix[col][sub_col] );                        
                        };
                    };   
                };      
            };
            return matrix;
        };

        vector<double> matrix_vector_multiply(vector<vector<double>> inputMatrix, vector<double> inputVector)
        {
            vector<double> product;
            for (int row = 0; row < (int)inputMatrix.size(); row++)
            {
                double sum = 0.0;
                for (int col = 0; col < (int)inputMatrix[0].size(); col++)
                {
                    sum += inputMatrix[row][col] * inputVector[col];
                }
                product.push_back(sum);
            }
            return product;
        };

        vector<double> get_direction(string method, 
                            bool condition, int kth_iteration, 
                            vector<vector<double>> Conditioner_factored, 
                            vector<double> residual, 
                            vector<double> r_prior, 
                            vector<double> p_prior
                            )
        {
            vector<double> p = p_prior;
            if (method == "RF") // Richardson Family
            {
                int len = (int)p.size();
                int index = kth_iteration % len;
                transform(p.begin(), p.end(), p.begin(), [](double i){return i=0;});
                p[index] = 1;
            }else if (method == "SD") // Steepest Decent
            {
                p = residual;
            }else if (method == "CG") // Conjugate Gradient
            {
                if (kth_iteration == 0)
                {
                    p = residual;
                }else
                {
                    double gamma = inner_product(residual,residual) / inner_product(r_prior, r_prior);
                    transform(p.begin(), p.end(), p.begin(), [gamma](double &c){ return c*gamma; }); // p = gamma * p
                    transform(p.begin(), p.end(), residual.begin(), p.begin(), plus<double>()); // p = r + gamma*p
                }
                
            };   

            if (condition) // Apply PreConditioner
            {
                p = solveCholesky(Conditioner_factored, p);
            };
            return p;
        };

        double inner_product(vector<double> first, vector<double> second)
        {
            double val = 0.0;
            for (int i = 0; i < (int)first.size(); i++)
            {
                val += first[i] * second[i];
            }
            return val;
        };

        vector<double> solveCholesky(vector<vector<double>> conditioner, vector<double> p)
        {
            int row = (int)conditioner.size();
            // Solve Ly=p
            vector<double> Solution;
            Solution.resize(p.size());
            for (int i = 0; i < row; i++)
            {
                Solution[i] = p[i];
                for (int j = i-1; j >= 0; j--)
                {
                    
                    Solution[i] -= Solution[j]*conditioner[i][j];
                    
                }
                Solution[i] /= conditioner[i][i];
            };
            // Solve L^Tp=y
            for (int i = row-1; i >= 0; i--)
            {
                for (int j = i+1; j < row; j++)
                {
                    // Because Upper Right is Transpose Lower Left we switch ith row with jth column for Lower val
                    Solution[i] -= Solution[j]*conditioner[j][i];
                    
                }
                Solution[i] /= conditioner[i][i];
            } 
            return Solution;     
        };

        void solve_system(string method, bool precondition)
        {
            includePrecondition = precondition;
            solve_exact_via_richardson_gauss_seidel();
            if (method == "RF")
            {
                // solve_exact_via_richardson_gauss_seidel();
                solve_via_richardson_gauss_seidel();
            }else if (method == "SD")
            {
                // solve_exact_via_steepest_descent();
                solve_via_steepest_descent();
            }else if (method == "CG")
            {
               //  solve_exact_via_conjugate_gradient();
                solve_via_conjugate_gradient();
            }else if (method == "GS")
            {
                solve_via_gauss_southwell();
            }else if (method == "SD_slow")
            {
                solve_via_slow_steepest_descent();
            } else if (method == "CG_newton")
            {
                solve_via_newton_trunc_CG();
            } else if (method == "BGFS")
            {
                solve_via_quasi_newt_BFGS();
            } 
        };

        void print_matrix(vector<vector<float>> matrix)
        {
            int row = (int)matrix.size();
            int column = (int)matrix[0].size();

            for(int i = 0; i < row; i++)
            {
                for(int j = 0; j < column; j++)
                {
                    cout << matrix[i][j] << " ";
                    if (j == column-1)
                    {
                        cout << "\n";
                    };
                }
            }
            
        };

        void print_solution(string method)
        {
            if (method == "RF")
            {
                cout << "RF Solution Vector \n";
                for (int i=0; i<(int)solution_vector_sp_RF.size(); i++)
                {
                    cout << solution_vector_sp_RF[i] << "\n";
                };

            }else if (method == "SD")
            {
                cout << "SD Solution Vector \n";
                for (int i=0; i<(int)solution_vector_sp_SD.size(); i++)
                {
                    cout << solution_vector_sp_SD[i] << "\n";
                };
            }else if (method == "SD_slow")
            {
                cout << "SD_slow Solution Vector \n";
                for (int i=0; i<(int)solution_vector_sp_SD_slow.size(); i++)
                {
                    cout << solution_vector_sp_SD_slow[i] << "\n";
                };
            }else if (method == "GS")
            {
                cout << "GS Solution Vector \n";
                for (int i=0; i<(int)solution_vector_sp_GS.size(); i++)
                {
                    cout << solution_vector_sp_GS[i] << "\n";
                };
            }else if (method == "CG")
            {
                cout << "CG Solution Vector \n";
                for (int i=0; i<(int)solution_vector_sp_CG.size(); i++)
                {
                    cout << solution_vector_sp_CG[i] << "\n";
                };
            }else if (method == "CG_newton")
            {
                cout << "Newton Truncated CG Solution Vector \n";
                for (int i=0; i<(int)solution_vector_sp_CG_newton.size(); i++)
                {
                    cout << solution_vector_sp_CG_newton[i] << "\n";
                };
            }else if (method == "BGFS")
            {
                cout << "Quasi Newton BGFS Solution Vector \n";
                for (int i=0; i<(int)solution_vector_sp_quasi_newt_BFGS.size(); i++)
                {
                    cout << solution_vector_sp_quasi_newt_BFGS[i] << "\n";
                };
            }
            
            
        };

        // Single Precision
        vector<float> Iterate(
            string method, 
            vector<vector<float>> A_matrix, 
            vector<float> b_vector, 
            vector<float> x_vector, 
            vector<vector<float>> preConditioner,
            bool condition
            )
        {
            vector<float> vector_sol = x_vector;
            vector<vector<float>> Conditioner_factored = cholesky(preConditioner);
            vector<float> conditioned_bVector = solveCholesky(Conditioner_factored, b_vector);
            vector<float> dummy_vector = matrix_vector_multiply(A_matrix, x_vector);
            vector<float> residual = b_vector;
            vector<float> r_prior = residual;
            vector<float> p = residual;
            transform(residual.begin(), residual.end(), dummy_vector.begin(), residual.begin(), minus<float>());
            double norm_two_b = norm_two(b_vector);
            double norm_two_conditioned_b = norm_two(conditioned_bVector);
            double norm_two_exact = norm_two_dp(solution_vector_dp);
            double metric;
            int k = 0;

            if (method != "CG_newton" || method != "BGFS")
            {
                bool zeros = all_of(residual.begin(), residual.end(), [](float i) { return i==0.0; });
                while (zeros == false && k<=30)
                {
                    p = get_direction(method, condition, k, Conditioner_factored, residual, r_prior, p);
                    vector<float> q = matrix_vector_multiply(A_matrix, p);
                    float alpha = inner_product(p, residual) / inner_product(p, q);
                    if (method=="SD_slow")
                    {
                        alpha = 1.1*alpha;
                    };
                    transform(p.begin(), p.end(), p.begin(), [alpha](float &c){ return c*alpha; }); // p = alpha * p
                    transform(vector_sol.begin(), vector_sol.end(), p.begin(), vector_sol.begin(), plus<float>()); // x = x + alpha*p
                    transform(q.begin(), q.end(), q.begin(), [alpha](float &c){ return c*alpha; }); // q = alpha * q
                    transform(residual.begin(), residual.end(), q.begin(), residual.begin(), minus<float>()); // r = r - alpha*q
                    // Calc and Store Kth Metrics -------------------------------
                    metric = norm_two_dp(vector_diff(vector_sol, solution_vector_dp)) / norm_two_exact;
                    metric_known.push_back(metric);
                    if (condition)
                    {
                        metric = norm_two(residual) / norm_two_conditioned_b;
                        metric_precondition.push_back(metric);
                    }else
                    {
                        metric = norm_two(residual) / norm_two_b;
                        metric_unknown.push_back(metric);
                    }
                    
                    // ------------------------------------------------------------ 
                    k++;
                    zeros = all_of(residual.begin(), residual.end(), [](float i) { return i==0.0; });
                };
            }else if (method == "CG_newton")
            {
                bool zeros = all_of(residual.begin(), residual.end(), [](float i) { return i==0.0; });
                while  (zeros == false && k<=50)
                {
                
                    // Truncated Conjugate Gradient Inner Loop
                    bool terminate = false;
                    vector<float> g = residual;

                    vector<float> d = g;
                    float neg = -1;
                    transform(d.begin(), d.end(), d.begin(), [neg](float &c){ return c*neg; });
                    
                    int size = vector_sol.size();
                    vector<float> z(size);
                    fill(z.begin(), z.end(), 0);

                    float epsilon = min(0.5, sqrt(norm_two(g)))*norm_two(g);

                    int i = 0; 
                    while (!terminate) 
                    {
                        vector<float> g_prior = g;
                        vector<float> q = matrix_vector_multiply(A_matrix, d);
                        float mu = inner_product(d, q);
                        if (mu <= 0)
                        {
                            if (i=0)
                            {
                                p = d;
                                break;
                            }else
                            {
                                p = z;
                                break;
                            };
                        };
                        float alpha = inner_product(g, g) / mu;
                        transform(d.begin(), d.end(), d.begin(), [alpha](float &c){ return c*alpha; }); // d = alpha * d
                        transform(z.begin(), z.end(), d.begin(), z.begin(), plus<float>()); // z = z + alpha*d
                        transform(q.begin(), q.end(), q.begin(), [alpha](float &c){ return c*alpha; }); // q = alpha * q
                        transform(g.begin(), g.end(), q.begin(), g.begin(), minus<float>()); // r = r - alpha*q
                        if (norm_two(g) <= epsilon)
                        {
                            p = z;
                            break;
                        };
                        
                        float beta = inner_product(g, g) / inner_product(g_prior, g_prior);
                        alpha = beta / alpha;
                        transform(d.begin(), d.end(), d.begin(), [alpha](float &c){ return c*alpha; }); // d = beta * d 
                        transform(d.begin(), d.end(), g.begin(), d.begin(), minus<float>()); // d = beta*d - g
                        i++;
                    };

                    // float sig = inner_product(residual, residual);
                    // vector<float> que = matrix_vector_multiply(A_matrix, p);
                    // float mew = inner_product(p, que);
                    // float alph = sig / mew;
                    // transform(p.begin(), p.end(), p.begin(), [alph](float &c){ return c*alph; }); // d = alpha * d
                    transform(vector_sol.begin(), vector_sol.end(), p.begin(), vector_sol.begin(), plus<float>()); // x = x + alpha*p
                    // transform(residual.begin(), residual.end(), que.begin(), residual.begin(), minus<float>()); // r = r - alpha*q

                    // Calc and Store Kth Metrics -------------------------------
                    metric = norm_two_dp(vector_diff(vector_sol, solution_vector_dp)) / norm_two_exact;
                    metric_known.push_back(metric);
                    if (condition)
                    {
                        metric = norm_two(residual) / norm_two_conditioned_b;
                        metric_precondition.push_back(metric);
                    }else
                    {
                        metric = norm_two(residual) / norm_two_b;
                        metric_unknown.push_back(metric);
                    }
                    k++;
                    zeros = all_of(residual.begin(), residual.end(), [](float i) { return i==0.0; });
                };

            }else if (method == "BGFS")
            {
                vector<float> s;
                while (sqrt(norm_two(s)) <= 0.001 && k<=50)
                {
                    p = get_direction(method, condition, k, Conditioner_factored, residual, r_prior, p); // solve B_k * p_k = -f'(x_k) via CG
                    vector<float> q = matrix_vector_multiply(A_matrix, p);
                    float alpha = inner_product(p, residual) / inner_product(p, q);
                    vector<float> y = b_vector;
                    transform(p.begin(), p.end(), p.begin(), [alpha](float &c){ return c*alpha; }); // d = alpha * d
                    s = p; // s = x_k+1 - x_k
                    vector<float> y = p; // y = f'(x_k+1) - f'(x_k)
                    transform(vector_sol.begin(), vector_sol.end(), p.begin(), vector_sol.begin(), plus<float>()); // x = x + alpha*p
                    Bmatrix = update_matrix(s, y, Bmatrix, "BFGS");
                    
                    k++;
                };

            };
            
            return vector_sol;
        };

        vector<vector<float>> cholesky(vector<vector<float>> matrix)
        {
            // Factor Matrix In Place
            int m_rows = (int)matrix.size();
            int n_columns = (int)matrix[0].size();

            for (int col = 0; col < n_columns; col++)
            {
                if (matrix[col][col] == 0)
                {
                    throw runtime_error("Math error: Attempted to divide by Zero. Consider Pivoting Matrix\n");
                }else
                {
                    // Root of Top Left
                    matrix[col][col] = sqrt(matrix[col][col]);
                    // Calculate Lower Trapeziodal Elements
                    for (int row = col+1; row < m_rows; row++)
                    {
                        matrix[row][col] = matrix[row][col] / matrix[col][col];
                    }
                    // Update SubMatrix Elements Before Next Iteration
                    for (int sub_row = col+1; sub_row < m_rows; sub_row++)
                    {
                        for (int sub_col = col+1; sub_col < n_columns; sub_col++)
                        {
                            matrix[sub_row][sub_col] = matrix[sub_row][sub_col] - ( matrix[sub_row][col] * matrix[col][sub_col] );                        
                        };
                    };   
                };      
            };
            return matrix;
        };

        vector<float> matrix_vector_multiply(vector<vector<float>> inputMatrix, vector<float> inputVector)
        {
            vector<float> product;
            for (int row = 0; row < (int)inputMatrix.size(); row++)
            {
                float sum = 0.0;
                for (int col = 0; col < (int)inputMatrix[0].size(); col++)
                {
                    sum += inputMatrix[row][col] * inputVector[col];
                }
                product.push_back(sum);
            }
            return product;
        };

        vector<float> get_direction(string method, 
                            bool condition, int kth_iteration, 
                            vector<vector<float>> Conditioner_factored, 
                            vector<float> residual, 
                            vector<float> r_prior, 
                            vector<float> p_prior
                            )
        {
            vector<float> p = p_prior;
            if (method == "RF") // Richardson Family
            {
                int len = (int)p.size();
                int index = kth_iteration % len;
                transform(p.begin(), p.end(), p.begin(), [](float i){ return i=0;});
                p[index] = 1;
            }else if (method == "SD") // Steepest Decent
            {
                p = residual;
            }else if (method == "SD_slow") // SLOW Steepest Decent
            {
                p = residual;
            }else if (method == "GS") // Gauss Southwell
            {
                p = residual;
                int i = index_of_max_mag(p);
                int sign = (p[i] > 0) - (p[i] < 0);
                transform(p.begin(), p.end(), p.begin(), [](float i){ return i=0;});
                p[i] = sign;

            }else if (method == "CG") // Conjugate Gradient
            {
                if (kth_iteration == 0)
                {
                    p = residual;
                }else
                {
                    float gamma = inner_product(residual,residual) / inner_product(r_prior, r_prior);
                    transform(p.begin(), p.end(), p.begin(), [gamma](float &c){ return c*gamma; }); // p = gamma * p
                    transform(p.begin(), p.end(), residual.begin(), p.begin(), plus<float>()); // p = r + gamma*p
                }
                
            };   

            if (condition) // Apply PreConditioner
            {
                p = solveCholesky(Conditioner_factored, p);
            };
            return p;
        };

        vector<vector<float>> update_matrix(vector<float> s, vector<float> y, vector<vector<float>> Bmatrix, string method)
        {
            
            return Bmatrix;
        };

        float inner_product(vector<float> first, vector<float> second)
        {
            float a,b;
            float val = 0.0;
            for (int i = 0; i < (int)first.size(); i++)
            {
                a =  first[i];
                b = second[i];
                val += a * b;
            }
            return val;
        };

        float index_of_max_mag(vector<float> vec)
        {
            int index = 0;
            float val = 0.0;
            for (int i = 0; i < (int)vec.size(); i++)
            {
                if (abs(vec[i]) > val)
                {
                    index = i;
                    val = abs(vec[i]);
                };
                
            };
            return index;
        };

        vector<float> solveCholesky(vector<vector<float>> conditioner, vector<float> p)
        {
            int row = (int)conditioner.size();
            // Solve Ly=p
            vector<float> Solution;
            Solution.resize(p.size());
            for (int i = 0; i < row; i++)
            {
                Solution[i] = p[i];
                for (int j = i-1; j >= 0; j--)
                {
                    
                    Solution[i] -= Solution[j]*conditioner[i][j];
                    
                }
                Solution[i] /= conditioner[i][i];
            };
            // Solve L^Tp=y
            for (int i = row-1; i >= 0; i--)
            {
                for (int j = i+1; j < row; j++)
                {
                    // Because Upper Right is Transpose Lower Left we switch ith row with jth column for Lower val
                    Solution[i] -= Solution[j]*conditioner[j][i];
                    
                }
                Solution[i] /= conditioner[i][i];
            } 
            return Solution;     
        };

        // Metrics
        double norm_two(vector<float> vect)
        {
            int row = (int)vect.size();
            double maxVal = 0.0;

            for(int i = 0; i < row; i++)
            {
                maxVal = maxVal + pow(abs(vect[i]),2);
            };           
            
            return sqrt(maxVal);
        };

        double norm_two_dp(vector<double> vect)
        {
            int row = (int)vect.size();
            double maxVal = 0.0;

            for(int i = 0; i < row; i++)
            {
                maxVal = maxVal + pow(abs(vect[i]),2);
            };           
            
            return sqrt(maxVal);
        };

        vector<double> vector_diff(vector<float> vect1, vector<double> vect2)
        {
            vector<double> difference;
            transform(vect1.begin(), vect1.end(),back_insert_iterator(difference),[](float val) -> double{return (double)val;});
            transform(difference.begin(), difference.end(), vect2.begin(), difference.begin(), minus<double>());
            return difference;
        };

        vector<double> get_metric_known()
        {
            return metric_known;
        };

        vector<double> get_metric_unknown()
        {
            return metric_unknown;
        };

        vector<double> get_metric_preconditioned()
        {
            return metric_precondition;
        };
};

int main()
{
    bool applyPreConditioner = false;
    string method = "CG_newton";
    int condition = 0;
    int pc_matrix = 3;

    vector<vector<double>> inputMatrix = readfile_matrix("data_matrix.txt");
    vector<vector<double>> inputPreconditioner = readfile_matrix("data_preCondition.txt");
    vector<double> inputVector_b = readFile_vector("data_vector_b.txt");
    vector<double> inputVector_x = readFile_vector("data_vector_x.txt");

    IterativeMethods Iterate1(inputMatrix, inputVector_b, inputVector_x, inputPreconditioner, applyPreConditioner);
    Iterate1.solve_system(method, applyPreConditioner);
    Iterate1.print_solution(method);
    vector<double> metric_known, metric_unknown, metric_conditioned;
    metric_known = Iterate1.get_metric_known();
    metric_unknown = Iterate1.get_metric_unknown();
    metric_conditioned = Iterate1.get_metric_preconditioned();

    // // Write Data to CSV File
    // ofstream myfile;
    // myfile.open("analysis_data.csv", ofstream::app);
    // if (applyPreConditioner)
    // {
    //     for (int i = 0; i < (int)metric_known.size(); i++)
    //     {
    //         myfile << "\n";
    //         myfile << method << "," << i << "," << metric_known[i] << "," << metric_conditioned[i] << "," << pc_matrix << "," << condition << ",";
    //     }
    // }else
    // {
    //     for (int i = 0; i < (int)metric_known.size(); i++)
    //     {
    //         myfile << "\n";
    //         myfile << method << "," << i << "," << metric_known[i] << "," << metric_unknown[i] << "," << 0 << "," << condition << ",";
    //     }
    // };
    // myfile.close();
 }