#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
using namespace std;

vector<float> readFile_vector(string filename){
    vector<float> data;
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

class Matrix
{
    private:
        // attributes
        vector<vector<float>> A_matrix;
        vector<vector<float>> LU_matrix;
        vector<vector<float>> LU_product_matrix;
        vector<vector<float>> A_product;
        vector<vector<float>> QR_matrix;
        vector<vector<float>> QR_product_matrix;
        vector<float> tau;
        vector<float> Solution_LU;
        vector<float> Solution_QR;
        vector<float> b_vector;
        vector<float> b_vector_approx;
        vector<float> least_square_solution;
        vector<int> Permute_row;
        vector<int> Permute_column;

        // methods
        void LU_Gauss()
        {
            int m_rows = (int)LU_matrix.size();
            int n_columns = (int)LU_matrix[0].size();

            for (int col = 0; col < n_columns; col++)
            {
                if (LU_matrix[col][col] == 0)
                {
                    throw runtime_error("Math error: Attempted to divide by Zero. Consider Pivoting Matrix\n");
                }else
                {
                    // Calculate Lower Trapeziodal Elements
                    for (int row = col+1; row < m_rows; row++)
                    {
                        LU_matrix[row][col] = LU_matrix[row][col] / LU_matrix[col][col];
                    }
                    // Update SubMatrix Elements Before Next Iteration
                    for (int sub_row = col+1; sub_row < m_rows; sub_row++)
                    {
                        for (int sub_col = col+1; sub_col < n_columns; sub_col++)
                        {
                            LU_matrix[sub_row][sub_col] = LU_matrix[sub_row][sub_col] - ( LU_matrix[sub_row][col] * LU_matrix[col][sub_col] );                        
                        };
                    };   
                };      
            };
        };

        void LU_PartialPivot()
        {
            int m_rows = (int)LU_matrix.size();
            int n_columns = (int)LU_matrix[0].size();

            for (int col = 0; col < n_columns; col++)
            {
                partial_pivot(col, m_rows);
                if (LU_matrix[col][col] == 0)
                {
                    throw runtime_error("Math error: Attempted to divide by Zero. Consider Pivoting Matrix\n");
                }else
                {
                    // Calculate Lower Trapeziodal Elements
                    for (int row = col+1; row < m_rows; row++)
                    {
                        LU_matrix[row][col] = LU_matrix[row][col] / LU_matrix[col][col];
                    }
                    // Update SubMatrix Elements Before Next Iteration
                    for (int sub_row = col+1; sub_row < m_rows; sub_row++)
                    {
                        for (int sub_col = col+1; sub_col < n_columns; sub_col++)
                        {
                            LU_matrix[sub_row][sub_col] = LU_matrix[sub_row][sub_col] - ( LU_matrix[sub_row][col] * LU_matrix[col][sub_col] );                        
                        };
                    };   
                };      
            };

        };

        void LU_CompletePivot()
        {
            int m_rows = (int)LU_matrix.size();
            int n_columns = (int)LU_matrix[0].size();

            for (int col = 0; col < n_columns; col++)
            {
                vector<int> maxVal_coordinates = find_maxVal(col, col);
                swap_columns(col, maxVal_coordinates[1]);
                swap_rows(col, maxVal_coordinates[0]);

                if (LU_matrix[col][col] == 0)
                {
                    throw runtime_error("Math error: Attempted to divide by Zero. Consider Pivoting Matrix\n");
                }else
                {
                    // Calculate Lower Trapeziodal Elements
                    for (int row = col+1; row < m_rows; row++)
                    {
                        LU_matrix[row][col] = LU_matrix[row][col] / LU_matrix[col][col];
                    }
                    // Update SubMatrix Elements Before Next Iteration
                    for (int sub_row = col+1; sub_row < m_rows; sub_row++)
                    {
                        for (int sub_col = col+1; sub_col < n_columns; sub_col++)
                        {
                            LU_matrix[sub_row][sub_col] = LU_matrix[sub_row][sub_col] - ( LU_matrix[sub_row][col] * LU_matrix[col][sub_col] );                        
                        };
                    };   
                };      
            };
        };

        vector<int> find_maxVal(int start_row, int start_column)
        {
            int row = start_row, column = start_column;
            float maxval = LU_matrix[start_row][start_column];
            for (int i = start_row; i < (int)LU_matrix.size(); i++)
            {
                for (int j = start_column; j < (int)LU_matrix[0].size(); j++)
                {
                    if (maxval < LU_matrix[i][j])
                    {
                        maxval = LU_matrix[i][j];
                        row = i;
                        column = j;
                    }
                    
                }
                
            }
            
            vector<int> maxVal_loc = {row, column};
            return maxVal_loc;
        };

        void swap_columns(int initial_col, int maxval_col)
        {
            // Update Matrix
            if (initial_col != maxval_col)
            {
                float val;
                int num_rows = (int)LU_matrix[0].size();
                for (int i = 0; i < num_rows; i++)
                {
                    val = LU_matrix[i][maxval_col];
                    LU_matrix[i][maxval_col] = LU_matrix[i][initial_col];
                    LU_matrix[i][initial_col] = val;
                }
                
                // Update Permutation Vector
                int value;
                value = Permute_column[maxval_col];
                Permute_column[maxval_col] = Permute_column[initial_col];
                Permute_column[initial_col] = value;
            }
        };

        void swap_rows(int initial_row, int maxval_row)
        {
            // Update Matrix
            if (initial_row != maxval_row)
            {
                float val;
                int num_columns = (int)LU_matrix[0].size();
                for (int i = 0; i < num_columns; i++)
                {
                    val = LU_matrix[maxval_row][i];
                    LU_matrix[maxval_row][i] = LU_matrix[initial_row][i];
                    LU_matrix[initial_row][i] = val;
                }
                
                // Update Permutation Vector
                int value;
                value = Permute_row[maxval_row];
                Permute_row[maxval_row] = Permute_row[initial_row];
                Permute_row[initial_row] = value;
            }

        };

        void create_row_vector(int num_rows)
        {
            for (int i = 0; i < num_rows; i++)
            {
                Permute_row.push_back(i);
            }
        };

        void create_col_vector(int num_cols)
        {
            for (int i = 0; i < num_cols; i++)
            {
                Permute_column.push_back(i);
            }
        };

        void partial_pivot(int col_num, int num_rows)
        {
            // Find Row with max value
            vector<float> columnVals;
            for (int i = col_num; i < num_rows; i++)
            {
                columnVals.push_back(LU_matrix[i][col_num]);
            }

            int row = col_num;
            float maxval = columnVals[0];
            for (int i = 0; i < (int)columnVals.size(); i++)
            {
                if (maxval < columnVals[i])
                {
                    maxval = columnVals[i];
                    row = col_num + i;
                }
                
            }

            // Pivot if there is a max value else where
            if (row > col_num)
            {
                float row_val;
                for (int i = 0; i < (int)LU_matrix[0].size(); i++)
                {
                    row_val = LU_matrix[row][i];
                    LU_matrix[row][i] = LU_matrix[col_num][i];
                    LU_matrix[col_num][i] = row_val;
                }
            // Store Row Permuations
                int index;
                index = Permute_row[row];
                Permute_row[row] = Permute_row[col_num];
                Permute_row[col_num] = index;
            }
                
        };

        vector<vector<float>> LU_product()
        {
            int num_rows = (int)LU_matrix.size();
            int num_columns = (int)LU_matrix[0].size();
            float Uval, Lval, sum;
            vector<vector<float>> product(num_rows, vector<float>(num_columns));

            for (int row = 0; row < num_rows; row++)
            {
                for (int col = 0; col < num_columns; col++)
                {
                    sum = 0;
                    for (int i = 0; i < num_columns; i++)
                    {
                        Uval = get_Uval(i, col);
                        Lval = get_Lval(row, i);
                        sum = sum + ( Uval * Lval );
                    };
                    
                    product[row][col] = sum;
                }
                
            }
            

            return product;
        };

        float get_Uval(int row, int column)
        {
            float Uval;

            if (row > column)
            {
                Uval = 0.0;
            }else
            {
                Uval = LU_matrix[row][column];
            }
            
            return Uval;
        };

        float get_Lval(int row, int column)
        {
            float Lval;

            if (row < column)
            {
                Lval = 0.0;
            }else if (row == column)
            {
                Lval = 1.0;
            }else
            {
                Lval = LU_matrix[row][column];
            }
            
            return Lval;
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

        void QR_Householder()
        {
            int m_rows = (int)QR_matrix.size();
            int n_columns = (int)QR_matrix[0].size();

            for (int col = 0; col < n_columns; col++)
            {
                if (QR_matrix[col][col] == 0)
                {
                    throw runtime_error("Math error: Attempted to divide by Zero. Consider Pivoting Matrix\n");
                }else
                {
                    // Set Lower Trapeziodal Elements ( Householder Vectors ) and construct tau_vector
                    set_Householder_vector(col, m_rows);

                    vector<float> w = set_w_vector(col, col+1, n_columns);

                    set_UpperRightRow(col, col+1, n_columns, w);

                    // Update SubMatrix Elements Before Next Iteration
                    set_Submatrix(col+1, col+1, n_columns, w);   
                };      
            };
            // Copy QR Matrix to matrix "QR_product_matrix" that will be used to reconstruct A
            copy(QR_matrix.begin(), QR_matrix.end(),back_inserter(QR_product_matrix));
        };

        void set_Householder_vector(int column, int numRows)
        {
            // Get Column Vector
            vector<float> columnVector;
            for (int i = column; i < numRows; i++)
            {
                columnVector.push_back(QR_matrix[i][column]);
            };
            
            // Apply Reflector and Store Householder vector inplace of lower triangular zero elements
            float norm = norm_two(columnVector);
            float v;
            vector<float> u_vector;

            for (int i = column; i < numRows; i++)
            {
                if (i == column)
                {
                    float rho = -1 * (QR_matrix[i][column] / abs(QR_matrix[i][column])) * norm; // -sign(a_ii)*||x||_2
                    v = QR_matrix[i][column] - rho;
                    QR_matrix[i][column] =  rho;
                }else
                {
                    QR_matrix[i][column] = QR_matrix[i][column] / v;
                    u_vector.push_back(QR_matrix[i][column]);
                }
            };
            norm = norm_two(u_vector);
            float tauVal = (1 + pow(norm,2))/2; // (1 + uTu)/2
            tau.push_back(tauVal);
            
        };

        vector<float> set_w_vector(int startRow, int startColumn, int EndColumn)
        {
            // Creates (a_12T + u_21TA_22)/tau
            vector<float> w;
            for (int col = startColumn; col < EndColumn; col++)
            {
                float wVal = QR_matrix[startRow][col];
                for (int row = startColumn; row < EndColumn; row++)
                {
                    wVal += QR_matrix[row][startColumn - 1] * QR_matrix[row][col];
                }
                wVal /= tau[startColumn - 1];
                w.push_back(wVal);
            }
            
            return w;
        };

        void set_UpperRightRow(int startRow, int startColumn, int EndColumn, vector<float> w)
        {
            int j=0;
            for (int i = startColumn; i < EndColumn; i++)
            {
                QR_matrix[startRow][i] -= w[j];
                j++;
            }
            
        };

        void set_Submatrix(int startRow, int startColumn, int EndColumn, vector<float> w)
        {
            int i=0;
            for (int col = startColumn; col < EndColumn; col++)
            {
                for (int row = startRow; row < EndColumn; row++)
                {
                    QR_matrix[row][col] -= QR_matrix[row][startColumn - 1] * w[i];
                }
                i++;
            }
            
        };

        void QR_product()
        {
            int num_rows = (int)QR_matrix.size();
            int num_columns = (int)QR_matrix[0].size();
            vector<vector<float>> product(num_rows, vector<float>(num_columns));

            for (int col = num_columns - 1; col >= 0; col--)
            {
                apply_householder_reflector(col, num_columns);
            }
        };

        void apply_householder_reflector(int start_column, int max_column)
        {
            vector<float> w_vector;
            for (int col = start_column; col < max_column; col++)
            {
                float wVal = set_w(start_column, col, max_column) / tau[start_column];
                w_vector.push_back(wVal);
            }

            int i = (int)w_vector.size()-1;
            for (int col = max_column-1; col >= start_column; col--)
            {
                for (int row = max_column-1; row >= start_column; row--)
                {
                    if (row == start_column)
                    {
                        QR_product_matrix[row][col] -= w_vector[i];

                    }else if ((col == start_column) && (row != start_column))
                    {

                        QR_product_matrix[row][col] *= -w_vector[i];

                    }else
                    {
            
                        QR_product_matrix[row][col] -= QR_product_matrix[row][start_column] * w_vector[i];
                    }
                }
                i--;
            }
            
        };

        float set_w(int start_row, int current_column, int end_column)
        {
            float wVal = 0.0;
            if (start_row == current_column)
            {
                wVal = QR_product_matrix[start_row][current_column];
            }else
            {
                for (int row = start_row; row < end_column; row++)
                {  
                    
                    if (row == start_row)
                    {
                        wVal = QR_product_matrix[row][current_column];
                    }else
                    {
                        wVal += QR_product_matrix[row][start_row] * QR_product_matrix[row][current_column];
                    }   
                }
            };
            
            return wVal;
        };

        double norm_one(vector<vector<float>> matrix)
        {
            int row = (int)matrix.size();
            int column = (int)matrix[0].size();
            double maxVal = 0.0;

            for(int j = 0; j < column; j++)
            {
                double col_sum = 0.0;
                for(int i; i < row; i++)
                {
                    col_sum = col_sum + abs(matrix[i][j]);
                };

                if (col_sum >= maxVal)
                {
                    maxVal = col_sum;
                };
                
            };
            return maxVal;
        };

        double norm_frob(vector<vector<float>> matrix)
        {
            int row = (int)matrix.size();
            int column = (int)matrix[0].size();
            double maxVal = 0.0;

            for(int j = 0; j < column; j++)
            {
                for(int i = 0; i < row; i++)
                {
                    double ab = abs(matrix[i][j]);
                    maxVal = maxVal + pow(ab,2);
                };           
            };
            return sqrt(maxVal);
        };

        double norm_inf(vector<vector<float>> matrix)
        {
            int row = (int)matrix.size();
            int column = (int)matrix[0].size();
            double maxVal = 0.0;

            for(int i = 0; i < row; i++)
            {
                double row_sum = 0.0;
                for(int j; j < column; j++)
                {
                    row_sum = row_sum + abs(matrix[i][j]);
                };

                if (row_sum >= maxVal)
                {
                    maxVal = row_sum;
                };
                
            };
            return maxVal;
        };

        double norm_two(vector<float> vect)
        {
            int row = (int)vect.size();
            double maxVal = 0.0;

            for(int i; i < row; i++)
            {
                maxVal = maxVal + pow(abs(vect[i]),2);
            };           
            
            return sqrt(maxVal);
        };

        void solve_via_LU()
        {
            int row = (int)LU_matrix.size();
            // Solve Ly=b
            Solution_LU = b_vector;
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < i; j++)
                {
                    Solution_LU[i] -= Solution_LU[j]*LU_matrix[i][j];
                }
                
            };
            // Solve Ux=y
            for (int i = row-1; i >= 0; i--)
            {
                for (int j = i+1; j < row; j++)
                {
                    
                    Solution_LU[i] -= Solution_LU[j]*LU_matrix[i][j];
                    
                }
                Solution_LU[i] /= LU_matrix[i][i];
            }      
        };

        void solve_via_QR()
        {
            int row = (int)QR_matrix.size();
            int max_columns = (int)QR_matrix[0].size();
            Solution_QR = b_vector;
            // Solve Rx = QTb = H_n......H_1b
            for (int col = 0; col < max_columns; col++)
            {   // Calculate w = uTb / tau
                float wVal = Solution_QR[col];
                for (int i = col+1; i < row; i++)
                {
                    wVal += QR_matrix[i][col] * Solution_QR[i];
                }
                wVal /= tau[col];

                // Update b vector elements
                for (int i = col; i < row; i++)
                {
                    if (i == col)
                    {
                        Solution_QR[i] -= wVal;
                    }else
                    {
                        Solution_QR[i] -= QR_matrix[i][col] * wVal;
                    }
                }
                
            }

            // Rx = y
            for (int i = row-1; i >= 0; i--)
            {
                for (int j = i+1; j < row; j++)
                {
                    
                    Solution_QR[i] -= Solution_QR[j]*QR_matrix[i][j];
                    
                }
                Solution_QR[i] /= QR_matrix[i][i];
            } 
        };

        void compute_b(string factor_method)
        {
            int num_rows = A_matrix.size();
            int num_columns = A_matrix[0].size();
            vector<float> x_vector;

            // Get Solution Vector
            if (factor_method == "LU")
            {
                x_vector = Solution_LU;
            }else if (factor_method == "QR")
            {
                x_vector = Solution_QR;
            }
            
            // Compute Ax and to get b_approx
            for (int row = 0; row < num_rows; row++)
            {
                double bVal = 0.0;
                for (int col = 0; col < num_columns; col++)
                {
                    bVal += A_matrix[row][col] * x_vector[col];
                }
                b_vector_approx.push_back(bVal);
            }
        };

    public:     
        // Overloaded Constructors
        Matrix(vector<vector<float>> inputMatrix, vector<float> inputVector)
        {
            copy(inputMatrix.begin(), inputMatrix.end(),back_inserter(A_matrix));
            copy(inputMatrix.begin(), inputMatrix.end(),back_inserter(LU_matrix));
            copy(inputMatrix.begin(), inputMatrix.end(),back_inserter(QR_matrix));
            create_col_vector((int)LU_matrix.size());
            create_row_vector((int)LU_matrix[0].size());
            b_vector = inputVector;
        };

        // methods
        void LU_Factorization(string str_pivot = "none")
        {
            if(str_pivot == "none" )
            {
                try
                {
                    LU_Gauss();
                }
                catch(runtime_error& e)
                {
                    cout << e.what() << '\n';
                }
                
            }else if (str_pivot == "partial")
            {
                LU_PartialPivot();
            }else if (str_pivot == "complete")
            {
                LU_CompletePivot();
            }      
        };

        void QR_Decomposition()
        {
            QR_Householder();
        };

        void product(string productType)
        {
            if (productType == "LU")
            {
                LU_product_matrix  = LU_product();
            }else if (productType == "QR")
            {
                QR_product();
            }
            
        };

        void print_product_LU()
        {
            cout << "LU Product \n";
            print_matrix(LU_product_matrix);
        };

        void print_product_QR()
        {
            cout << "QR Product \n";
            print_matrix(QR_product_matrix);
        };

        void print_LU()
        {
            cout << "LU Matrix \n";
            print_matrix(LU_matrix);
        };

        void print_QR()
        {
            cout << "QR Matrix \n";
            print_matrix(QR_matrix);
        };

        void print_A()
        {
            cout << "A Matrix \n";
            print_matrix(A_matrix);
        };

        void print_Permutation_row()
        {
            cout << "Rows Permution Vector \n";
            for (int i=0; i<(int)Permute_row.size(); i++)
            {
                cout << Permute_row[i] << "\n";
            };
        };

        void print_Permutation_column()
        {
            cout << "Columns Permution Vector \n";
            for (int i=0; i<(int)Permute_column.size(); i++)
            {
                cout << Permute_column[i] << "\n";
            };
        };

        void solve_system(string method)
        {
            if (method == "LU")
            {
                solve_via_LU();

            }else if (method == "QR")
            {
                solve_via_QR();
            }
        };

        vector<double> accuracy_residual(string factor_method)
        {
            int length = (int)b_vector.size();
            double numerator = 0.0, denominator = 0.0, accVal = 0.0;
            vector<float> residual;
            vector<double> accuracy;

            // Compute appoximate b vector in Ax_tild = b_approx
            if (factor_method == "LU")
            {
                compute_b("LU"); 
            }else if (factor_method == "QR")
            {
                compute_b("QR");
            };

            // Calculate residual
            for (int i = 0; i < length; i++)
                {
                    residual.push_back(b_vector[i] - b_vector_approx[i]);
                } 

            // Measure accuracy accross multiple norms
            for (int i = 1; i < 4; i++)
            {
                switch (i)
                {
                    case 1:
                    {   
                        int length = b_vector.size();
                        vector<vector<float>> residual_matrix, b_matrix;
                        residual_matrix.resize(length, vector<float>(1));
                        b_matrix.resize(length, vector<float>(1));

                        for (int row = 0; row < length; row++)
                        {
                            residual_matrix[row][0] = residual[row];
                            b_matrix[row][0] = b_vector[row];
                        }
                        
                        numerator = norm_one(residual_matrix);
                        denominator = norm_one(b_matrix);
                        accVal = numerator / denominator;
                        break;
                    }
                    case 2:
                    {
                        numerator = norm_two(residual);
                        denominator = norm_two(b_vector);
                        accVal = numerator / denominator;
                        break;
                    }
                    case 3:
                    {
                        int length = b_vector.size();
                        vector<vector<float>> residual_matrix, b_matrix;
                        residual_matrix.resize(length, vector<float>(1));
                        b_matrix.resize(length, vector<float>(1));

                        for (int row = 0; row < length; row++)
                        {
                            residual_matrix[row][0] = residual[row];
                            b_matrix[row][0] = b_vector[row];
                        }

                        numerator = norm_inf(residual_matrix);
                        denominator = norm_inf(b_matrix);
                        accVal = numerator / denominator;
                        break;
                    }
                }

                accuracy.push_back(accVal);
            }

            return accuracy;
        };

        vector<double> accuracy_exact(string factor_method, vector<float> exact_solution)
        {
            int length = (int)exact_solution.size();
            double numerator = 0.0, denominator = 0.0, accVal = 0.0;
            vector<float> solution_computed, solution_exact, residual;
            vector<double> accuracy;

            // Compute appoximate b vector in Ax_tild = b_approx
            solution_exact = exact_solution;
            if (factor_method == "LU")
            {
                solution_computed = Solution_LU; 
            }else if (factor_method == "QR")
            {
                solution_computed = Solution_QR;
            };

            // Calculate residual
            for (int i = 0; i < length; i++)
                {
                    residual.push_back(solution_exact[i] - solution_computed[i]);
                } 

            // Measure accuracy accross multiple norms
            for (int i = 1; i < 4; i++)
            {
                switch (i)
                {
                    case 1:
                    {   
                        int length = exact_solution.size();
                        vector<vector<float>> residual_matrix, solution_exact_matrix;
                        residual_matrix.resize(length, vector<float>(1));
                        solution_exact_matrix.resize(length, vector<float>(1));

                        for (int row = 0; row < length; row++)
                        {
                            residual_matrix[row][0] = residual[row];
                            solution_exact_matrix[row][0] = solution_exact[row];
                        }
                        
                        numerator = norm_one(residual_matrix);
                        denominator = norm_one(solution_exact_matrix);
                        accVal = numerator / denominator;
                        break;
                    }
                    case 2:
                    {
                        numerator = norm_two(residual);
                        denominator = norm_two(solution_exact);
                        accVal = numerator / denominator;
                        break;
                    }
                    case 3:
                    {
                        int length = exact_solution.size();
                        vector<vector<float>> residual_matrix, solution_exact_matrix;
                        residual_matrix.resize(length, vector<float>(1));
                        solution_exact_matrix.resize(length, vector<float>(1));

                        for (int row = 0; row < length; row++)
                        {
                            residual_matrix[row][0] = residual[row];
                            solution_exact_matrix[row][0] = solution_exact[row];
                        }

                        numerator = norm_inf(residual_matrix);
                        denominator = norm_inf(solution_exact_matrix);
                        accVal = numerator / denominator;
                        break;
                    }
                }

                accuracy.push_back(accVal);
            }

            return accuracy;
        };

        vector<double> correctness(string factor_method)
        {
            vector<vector<float>> delta = A_matrix;
            int num_rows = (int)A_matrix.size();
            int num_columns = (int)A_matrix[0].size();
            vector<double> correct_ness;

            if (factor_method == "LU")
            {
                for (int row = 0; row < num_columns; row++)
                {
                    for (int col = 0; col < num_rows; col++)
                    {
                        delta[row][col] -= LU_product_matrix[row][col];
                    }
                    
                }
                
            }else if (factor_method == "QR")
            {
                for (int row = 0; row < num_columns; row++)
                {
                    for (int col = 0; col < num_rows; col++)
                    {
                        delta[row][col] -= QR_product_matrix[row][col];
                    }
                    
                }
            }

            double numerator = 0.0, denominator = 0.0, cVal = 0.0;
            for (int i = 1; i < 4; i++)
            {
                switch (i)
                {
                    case 1:
                    {
                        if (factor_method == "LU")
                        {
                            numerator = norm_one(delta);
                            denominator = norm_one(LU_product_matrix);
                            cVal = numerator / denominator;
                        }else if (factor_method == "QR")
                        {
                            numerator = norm_one(delta);
                            denominator = norm_one(QR_product_matrix);
                            cVal = numerator / denominator;
                        };
                        break;
                    }
                    case 2:
                    {
                        if (factor_method == "LU")
                        {
                            numerator = norm_frob(delta);
                            denominator = norm_frob(LU_product_matrix);
                            cVal = numerator / denominator;
                        }else if (factor_method == "QR")
                        {
                            numerator = norm_frob(delta);
                            denominator = norm_frob(LU_product_matrix);
                            cVal = numerator / denominator;
                        };
                        break;
                    }
                    case 3:
                    {
                        if (factor_method == "LU")
                        {
                            numerator = norm_inf(delta);
                            denominator = norm_inf(LU_product_matrix);
                            cVal = numerator / denominator;
                        }else if (factor_method == "QR")
                        {
                            numerator = norm_inf(delta);
                            denominator = norm_inf(LU_product_matrix);
                            cVal = numerator / denominator;
                        };
                        break;
                    }
                }

                correct_ness.push_back(cVal);
            }
            

            return correct_ness;
        };

        void solve_least_squares()
        {

        };

        void compute_growth_factor()
        {
            
        };

        void print_solution(string method)
        {
            if (method == "LU")
            {
                cout << "LU Solution Vector \n";
                for (int i=0; i<(int)Solution_LU.size(); i++)
                {
                    cout << Solution_LU[i] << "\n";
                };

            }else if (method == "QR")
            {
                cout << "QR Solution Vector \n";
                for (int i=0; i<(int)Solution_QR.size(); i++)
                {
                    cout << Solution_QR[i] << "\n";
                };
            }
        };

        vector<float> solution()
        {
            return Solution_LU;
        };

};


vector<float> generate_uniform_mesh(float start, float end, int numpoints)
{
    vector<float> mesh;
    float h = ( end - start ) / numpoints;
    // Generate points outside of mesh
    mesh.push_back(start - 3*h);
    mesh.push_back(start - 2*h);
    mesh.push_back(start-h);
    // Generate mesh
    for (int i = 0; i < numpoints; i++)
    {
        float val = start + (h * i);
        mesh.push_back(val);
    };
    mesh.push_back(end);
    // Generate points outside mesh
    mesh.push_back(end + h);
    mesh.push_back(end + 2*h);
    mesh.push_back(end + 3*h);

    return mesh;
};

void print_vector(vector<float> input_vector)
        {
            cout << "Vector Vals \n";
            for (int i=0; i<(int)input_vector.size(); i++)
            {
                cout << input_vector[i] << "\n";
            };
        };

float func1(float x)
{
    float y = pow(x,7) - pow(x,2);
    return y;
};

vector<float> slicing(vector<float>& arr,
                    int X, int Y)
{
 
    // Starting and Ending iterators
    auto start = arr.begin() + X;
    auto end = arr.begin() + Y + 1;
 
    // To store the sliced vector
    vector<float> result(Y - X + 1);
 
    // Copy vector using copy function()
    copy(start, end, result.begin());
 
    // Return the final sliced vector
    return result;
}

vector<float> get_observed_fvals(vector<float> mesh, int derivative)
{
    vector<float> observed_f;
    float del = pow(10,-9);

    if (derivative == 1)
    {
        float val = func1(mesh[0]+del) - func1(mesh[0]);
        val = val / del;
        observed_f.push_back(val);
    }else if (derivative == 2)
    {
        float val = func1(mesh[0]+del) - func1(mesh[0]-del);
        val = val / pow(del,2);
        observed_f.push_back(val);
    };

    for (int i = 0; i < (int)mesh.size(); i++)
    {
        float val = func1(mesh[i]);
        observed_f.push_back(val);
    }

    if (derivative == 1)
    {
        float val = func1(mesh.back()+del) - func1(mesh.back());
        val = val / del;
        observed_f.push_back(val);
    }else if (derivative == 2)
    {
        float val = func1(mesh.back()+del) - func1(mesh.back()-del);
        val = val / pow(del,2);
        observed_f.push_back(val);
    };

    return observed_f;
};

vector<vector<float>> Initialize_B_matrix(int dimension, float interval_width)
{
    vector<vector<float>> B_Matrix(dimension,vector<float>(dimension, 0.0));
    for (int row = 0; row < dimension; row++)
    {   
        if (row==0)
        {
            B_Matrix[row][row] = -3 / interval_width;
            B_Matrix[row][row+2] = 3 / interval_width;
        }else if (row==dimension-1)
        {
            B_Matrix[row][row-2] = -3 / interval_width;
            B_Matrix[row][row] = 3 / interval_width;
        }else
        {
            B_Matrix[row][row-1] = 1;
            B_Matrix[row][row] = 4;
            B_Matrix[row][row+1] = 1;
        };
    };
    return B_Matrix;
};

float cubic_piecewise_func(float x_fixed, float x_variable)
{
    float val=0.0;
    if (x_variable < x_fixed)
    {
        val = x_fixed-x_variable;
        val = pow(val,3);
    }
    return val;
};

float calculate_piecewise_cubic_polynomial(float xval, int index, vector<float> mesh)
{
    float val = cubic_piecewise_func(mesh[index+2], xval);
    val += -4 * cubic_piecewise_func(mesh[index+1], xval);
    val += 6 * cubic_piecewise_func(mesh[index], xval);
    val += -4 * cubic_piecewise_func(mesh[index-1], xval);
    val += cubic_piecewise_func(mesh[index-2], xval);
    return val;
};

vector<float> cubic_bspline_interpolation(vector<float> x_vector, vector<float> mesh, vector<float> alphas, float h)
{
    vector<float> interpolated_fvals;
    float norm_const = 1 / (pow(h,3)); // 1/ (3! * h^3 )

    for (int i = 0; i < (int)x_vector.size(); i++)
    {
        float fval = 0.0;
        for (int j = 2; j < (int)mesh.size()-2; j++)
        // starting at index 2 because that is where mesh truly starts prior 3 and last 3 are xtra points needed for bspline interpolation
        {
            float B_of_x = norm_const * calculate_piecewise_cubic_polynomial(x_vector[i], j, mesh);
            fval += alphas[j-2] * B_of_x;
        };
        interpolated_fvals.push_back(fval);
    };
    return interpolated_fvals;
};

int main()
{
    vector<float> x_to_interpolate  = readFile_vector("C:\\Users\\JustinCata.000\\Desktop\\FCM2_Projects\\Project7\\X_Vals.txt");

    float a = 0, b = 20, n=20; //a is start, b is end, n is number of points
    float h = ( b-a ) / n;
    vector<float> mesh = generate_uniform_mesh(a, b, n);
    vector<float> observed_f = get_observed_fvals(slicing(mesh,3,(int)mesh.size()-4), 1);
    vector<vector<float>> B_Matrix = Initialize_B_matrix((int)observed_f.size(), h);
    Matrix Matrix1(B_Matrix, observed_f);
    Matrix1.LU_Factorization();
    Matrix1.solve_system("LU");
    vector<float> alphas = Matrix1.solution();
    vector<float> f_interpolated = cubic_bspline_interpolation(x_to_interpolate, mesh, alphas, h);
    print_vector(f_interpolated);
};