#include "search_result.h"
#include "common.h"

std::ostream& operator<<(std::ostream& stream, const search_result& statistic) {
    return stream << "Search result = {Search method: " << (search_method_string.begin())[statistic.type] << "; " 
    << "Root of function: " << statistic.result << "; " 
    << "Accuracy: " << statistic.accuracy << "; " 
    << "Iterations: " << statistic.iterations << "; " 
    << "Function probes: " << statistic.function_probes << "}\n";
}