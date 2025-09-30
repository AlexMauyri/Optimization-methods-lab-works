#include "search_result_nd.h"
#include "common.h"
#include "custom.h"

std::ostream& operator<<(std::ostream& stream, const search_result_nd& statistic) {
    stream << "Search result = {Search method: " << (search_method_string_nd.begin())[statistic.type] << "; " 
    << "Extremum of function: "; 
    custom_vector_print(std::cout, statistic.result);
    stream << "; " 
    << "Accuracy: " << statistic.accuracy << "; " 
    << "Iterations: " << statistic.iterations << "; " 
    << "Function probes: " << statistic.function_probes << "}\n";

    return stream;
}