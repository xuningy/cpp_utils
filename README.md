# C++ Utils
#### Author: Xuning Yang <xuning@cmu.edu>
A set of C++ utility functions that I have written for research and development purposes.

-------------------
### `linalg_utils.h`: Defines a collection of linear algebra utility functions.
Given input vectors with type `T`:

Linspace: generates `N` linearly spaced values between lower bound (lb) and upper bound (ub), inclusive. If `lb==ub` and `N=1`, it returns a vector with 1 element (of that value).
```
std::vector<float> linearly_spaced_vector = Linspace(lb, ub, N);
std::vector<float> vector_of_one_element = Linspace(0.5, 0.5, 1);
```

-------------------
### `print_utils.h`: Defines a collection of print utility functions.

Prints a bunch of things: any standard sequential container (`std::vector`, `std::deque`, etc.), Eigen matrices wrt a specific layout.

```
std::vector<float> vec{10, 11, 12};
print(vec)                              // Output: (3): 10 11 12
print(vec, "name of vector")            // Output: name of vector (3): 10 11 12

Eigen::MatrixXd m;
m << 1, 2, 3,
     4, 5, 6,
     7, 8, 9;
print(m, "matrix")            // Output: matrix (3x3): [1, 2, 3;
                              //                        4, 5, 6;
                              //                        7, 8, 9]
```

-------------------
### `stats_utils.h`: Defines a collection of statistics utility functions.

RangeSample: draws `k` samples of  `T` uniformly from the range [lb, ub), and returns a vector. `T` is one of `{ int, float, double }`.

```
Eigen::VectorXd vec = RangeSample(lb, ub, k);
std::vector<float> vec = RangeSample(lb, ub, k);
```

DataSample: draws `k` samples sampled uniformly at random from `Container`
`data`, with replacement, and returns the sampled values in the same type
`Container`.

```
std::deque<double> sampled_vec = DataSample(data, k);
```

UniformDiscreteSample: draws `k` samples sampled from a discrete range from 0
to N, without replacement, and returns a `std::vector<int>`.

```
std::vector<int> sampled_numbers = UniformDiscreteSample(N, k);
```

DiscreteSample: `k` samples sampled at random according to distribution
`prob` with replacement, where `prob` is a probability array whose elements
sum to 1. T is one of { float, double }.

```
std::vector<int> sampled_numbers = UniformDiscreteSample(prob, k);
```

GaussianPdf: computes the probability density for a sample `x`, according to
mean `mean` and covariance `sigma`.

```
size_t M = 10;
Eigen::VectorXd mean(M);
Eigen::VectorXd x(M);
Eigen::MatrixXd Sigma(M,M);
float prob = GaussianPdf(x, mean, sigma);
```

Covariance: computes the sample covariance for an Eigen Matrix `samples` of
size MxN, and outputs an Eigen matrix of size MxM.
`samples` contains N samples of size M.

```
// computes covariance
Eigen::MatrixXd samples(M, N);
Eigen::MatrixXd sigma = Covariance(samples);              

// computes covariance about a given mean
Eigen::VectorXd mean(M);
Eigen::MatrixXd sigma = Covariance(samples, mean)

// computes weighted covariance about a given mean
Eigen::VectorXd weights(M);
Eigen::MatrixXd sigma = Covariance(samples, mean, weights)

```


-------------------
### `vector_utils.h`: Defines a collection of vector utility functions.
Given input vectors with type `T`:

```
std::vector<T> vec;
std::deque<T> deq;
```

Sort:
```
std::vector<T> vec_sorted;              // sorted vector will be stored here
std::vector<T> idx_sorted;              // index of the sorted vector will be stored here
Sort(vec, &vec_sorted, &idx_sorted)     
```
`vec` & `vec_sorted` must be of the same container type;  `idx_sorted` can be of any container type

Find:
```
T elem;
size_t index_of_found_element;   
bool element_found = Find(vec, elem);    
bool element_found = Find(vec, elem, &index_of_found_element);     
```

Max:
```
size_t index_of_max_element;
T max = Max(vec);
T max = Max(vec, &index_of_max_element);
```

Min:
```
size_t index_of_min_element;
T min = Min(vec);
T min = Min(vec, &index_of_min_element);
```

Invert: `vec` and `inverted_vec` can be of different container types
```
std::deque<T> inverted_vec;
Invert(vec, &inverted_vec);
```

Normalize: `vec` and `normalized_vec` can be of different container types
```
std::deque<T> normalized_vec;
Normalize(vec, &normalized_vec);
```
