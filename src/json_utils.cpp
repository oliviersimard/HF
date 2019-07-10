#include "json_utils.h"


const std::streampos Json_utils::getSize(const std::string& filename){

    std::streampos begin, end;
    std::ifstream myfile(filename, std::ifstream::in);
    begin = myfile.tellg();
    myfile.seekg(0, std::ifstream::end);
    end = myfile.tellg();
    myfile.close();

    return end-begin;
}

const json_spirit::mValue& Json_utils::get_object_item(const json_spirit::mValue& element, const std::string& name){

    return element.get_obj().at(name);
}

const json_spirit::mValue& Json_utils::get_array_item(const json_spirit::mValue& element, size_t index){

    return element.get_array().at(index);
}

MembCarrier::MembCarrier(double* dp, double* dp2, int* ip, bool* bp) : db_ptr(dp), db_ptr2(dp2), int_prt(ip), bool_ptr(bp){

}

MembCarrier Json_utils::JSONLoading(std::string& filename){

    std::ifstream JSONText(filename, std::ifstream::in);
    std::stringstream buffer;
    buffer << JSONText.rdbuf();
    json_spirit::mValue value;
    json_spirit::read(buffer,value); // Setting mValue object from buffer content.
    const std::string fileText(buffer.str()); // needs c_str() if printed using stdout.
    JSONText.close();
    std::cout << fileText.c_str() << std::endl;

    // Reading from json file
    //double precision
    const auto& U_val = get_object_item(value, "U");
    const auto& V_val = get_object_item(value, "V");
    const auto& beta_val = get_object_item(value, "beta");
    const auto& q_1D_val = get_object_item(value, "q_1D");
    //integers
    const auto& N_tau_val = get_object_item(value, "N_tau");
    const auto& N_it_val = get_object_item(value, "N_it");
    const auto& dims_val = get_object_item(value, "dims");
    const auto& gridK_val = get_object_item(value, "gridK");
    //bool
    const auto& precomK_val = get_object_item(value, "precomK");
    const auto& twoDMap_val = get_object_item(value, "2DMap");
    //array
    const auto& array_val = get_object_item(value, "q_2D");
    double container[2]; // Way to extract array inputs in json file!
    for (size_t i=0; i < array_val.get_array().size(); i++){
        container[i] = array_val.get_array().at(i).get_real();
    }

    std::cout << "Size of input file: " << getSize("params.json") << "\n";

    //Collecting the json variables to instantiate param object
    double U = U_val.get_real(); double V = V_val.get_real();
    double beta = beta_val.get_real(); double temperature = 1.0/beta;
    double q_1D = q_1D_val.get_real();
    int N_it = N_it_val.get_int(); int N_tau = N_tau_val.get_int();
    double d_tau = beta/N_tau; int dims = dims_val.get_int();
    int gridK = gridK_val.get_int();
    bool precomK = precomK_val.get_bool(); bool twoDMap = twoDMap_val.get_bool();

    //Creating the arrays
    double dub[6] = { U, V, beta, temperature, d_tau, q_1D };
    int integ[4] = { N_tau, N_it, dims, gridK };
    bool barr[2] = { precomK, twoDMap };

    return MembCarrier(dub,container,integ,barr);
}