# C++ Utils
#### Author: Xuning Yang <xuning@cmu.edu>
A set of C++ utility functions that I have written for research and development purposes.

## Examples
For the following examples, suppose:

```
std::vector<T> vec;
std::deque<T> deq;
```
and `T` can be any type.
###### `vector_utils.h`: Defines a collection of vector utility functions.
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
bool element_found = Find(vec, elem)    
bool element_found = Find(vec, elem, &index_of_found_element)     
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
