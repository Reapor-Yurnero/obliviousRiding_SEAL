// main.cpp
#include "seal/seal.h"
#include <algorithm>
#include <chrono>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

using namespace std;
using namespace seal;

/*
Helper function: Get poly_modulus_degree from a SEALContext.
*/
inline size_t get_N_from_context(const seal::SEALContext &context) {
    auto &context_data = *context.key_context_data();
    return context_data.parms().poly_modulus_degree();
}

/*
Helper function: Get plaintext modulus from a SEALContext.
*/
inline uint64_t get_T_from_context(const seal::SEALContext &context) {
    auto &context_data = *context.key_context_data();
    return context_data.parms().plain_modulus().value();
}

/*
Helper function: Prints the parameters in a SEALContext.
*/
inline void print_parameters(const seal::SEALContext &context)
{
    auto &context_data = *context.key_context_data();

    /*
    Which scheme are we using?
    */
    std::string scheme_name;
    switch (context_data.parms().scheme())
    {
    case seal::scheme_type::bfv:
        scheme_name = "BFV";
        break;
    case seal::scheme_type::ckks:
        scheme_name = "CKKS";
        break;
    default:
        throw std::invalid_argument("unsupported scheme");
    }
    std::cout << "/" << std::endl;
    std::cout << "| Encryption parameters :" << std::endl;
    std::cout << "|   scheme: " << scheme_name << std::endl;
    std::cout << "|   poly_modulus_degree: " << context_data.parms().poly_modulus_degree() << std::endl;

    /*
    Print the size of the true (product) coefficient modulus.
    */
    std::cout << "|   coeff_modulus size: ";
    std::cout << context_data.total_coeff_modulus_bit_count() << " (";
    auto coeff_modulus = context_data.parms().coeff_modulus();
    std::size_t coeff_modulus_size = coeff_modulus.size();
    for (std::size_t i = 0; i < coeff_modulus_size - 1; i++)
    {
        std::cout << coeff_modulus[i].bit_count() << " + ";
    }
    std::cout << coeff_modulus.back().bit_count();
    std::cout << ") bits" << std::endl;

    /*
    For the BFV scheme print the plain_modulus parameter.
    */
    if (context_data.parms().scheme() == seal::scheme_type::bfv)
    {
        std::cout << "|   plain_modulus: " << context_data.parms().plain_modulus().value() << std::endl;
    }

    std::cout << "\\" << std::endl;
}

/*
Helper function: Prints the `parms_id' to std::ostream.
*/
inline std::ostream &operator<<(std::ostream &stream, seal::parms_id_type parms_id)
{
    /*
    Save the formatting information for std::cout.
    */
    std::ios old_fmt(nullptr);
    old_fmt.copyfmt(std::cout);

    stream << std::hex << std::setfill('0') << std::setw(16) << parms_id[0] << " " << std::setw(16) << parms_id[1]
           << " " << std::setw(16) << parms_id[2] << " " << std::setw(16) << parms_id[3] << " ";

    /*
    Restore the old std::cout formatting.
    */
    std::cout.copyfmt(old_fmt);

    return stream;
}

/*
Helper function: Prints a vector of floating-point values.
*/
template <typename T>
inline void print_vector(std::vector<T> vec, std::size_t print_size = 4, int prec = 3)
{
    /*
    Save the formatting information for std::cout.
    */
    std::ios old_fmt(nullptr);
    old_fmt.copyfmt(std::cout);

    std::size_t slot_count = vec.size();

    std::cout << std::fixed << std::setprecision(prec);
    std::cout << std::endl;
    if (slot_count <= 2 * print_size)
    {
        std::cout << "    [";
        for (std::size_t i = 0; i < slot_count; i++)
        {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    else
    {
        vec.resize(std::max(vec.size(), 2 * print_size));
        std::cout << "    [";
        for (std::size_t i = 0; i < print_size; i++)
        {
            std::cout << " " << vec[i] << ",";
        }
        if (vec.size() > 2 * print_size)
        {
            std::cout << " ...,";
        }
        for (std::size_t i = slot_count - print_size; i < slot_count; i++)
        {
            std::cout << " " << vec[i] << ((i != slot_count - 1) ? "," : " ]\n");
        }
    }
    std::cout << std::endl;

    /*
    Restore the old std::cout formatting.
    */
    std::cout.copyfmt(old_fmt);
}

/*
Helper function: Prints a matrix of values.
*/
template <typename T>
inline void print_matrix(std::vector<T> matrix, std::size_t row_size)
{
    /*
    We're not going to print every column of the matrix (there are 2048). Instead
    print this many slots from beginning and end of the matrix.
    */
    std::size_t print_size = 5;

    std::cout << std::endl;
    std::cout << "    [";
    for (std::size_t i = 0; i < print_size; i++)
    {
        std::cout << std::setw(3) << std::right << matrix[i] << ",";
    }
    std::cout << std::setw(3) << " ...,";
    for (std::size_t i = row_size - print_size; i < row_size; i++)
    {
        std::cout << std::setw(3) << matrix[i] << ((i != row_size - 1) ? "," : " ]\n");
    }
    std::cout << "    [";
    for (std::size_t i = row_size; i < row_size + print_size; i++)
    {
        std::cout << std::setw(3) << matrix[i] << ",";
    }
    std::cout << std::setw(3) << " ...,";
    for (std::size_t i = 2 * row_size - print_size; i < 2 * row_size; i++)
    {
        std::cout << std::setw(3) << matrix[i] << ((i != 2 * row_size - 1) ? "," : " ]\n");
    }
    std::cout << std::endl;
}

/*
Helper function: Print line number.
*/
inline void print_line(int line_number)
{
    std::cout << "Line " << std::setw(3) << line_number << " --> ";
}

/*
Helper function: Compute euclidean distance between (a,b) and (c,d).
*/
inline uint64_t distance(uint64_t a, uint64_t b, uint64_t c, uint64_t d) {
    if (a > c) a, c = c, a;
    if (b > d) b, d = d, b;
    uint64_t x = c-a, y = d-b;
    // cout << a << " " << b << " " << c << " " << d << " " << x << " " << y << endl;
    return x*x + y*y; 
}


int main() {
    stringstream parms_stream; // Stream used to share parameters among server and clients.
    stringstream data_stream; // Stream used to share data (pk, rlk, ctxt) among server and clients.
    stringstream sk_stream; // Stream simulating the local secret storage of the rider.
    vector<vector<uint64_t> > driver_vectors; // Not used by anyone. For verification of correctness of the application.
    size_t driver_num; // Number of drivers in this region set by the server.

    /*
    Step 1: Server sets and sends the parameters.
    */
    {
        // Set the parameters.
        EncryptionParameters parms(scheme_type::bfv);
        size_t poly_modulus_degree = 8192;
        // size_t poly_modulus_degree = 4096;
        parms.set_poly_modulus_degree(poly_modulus_degree);
        parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
        uint64_t plain_modulus = 65929217;
        // uint64_t plain_modulus = 1032193;
        parms.set_plain_modulus(plain_modulus);

        // driver_num = 3;
        driver_num = poly_modulus_degree/2;

        auto size = parms.save(parms_stream);
        print_line(__LINE__);
        cout << "EncryptionParameters: wrote " << size << " bytes" << endl;
    }

    /*
    Step 2: The rider follows the parameters and generated a pair of pk and sk, 
    and then she uses either pk or sk to encrypt her coordinates to send to the server.
    */
    {
        // Load the parms from the parms_stream and reset the read pointer.
        EncryptionParameters parms;
        parms.load(parms_stream);
        parms_stream.seekg(0, parms_stream.beg);
        SEALContext context(parms);

        // Rider generetes key and stores them to stream properly.
        KeyGenerator keygen(context);
        SecretKey secret_key = keygen.secret_key();
        Serializable<PublicKey> public_key = keygen.create_public_key();
        Serializable<RelinKeys> relin_keys = keygen.create_relin_keys();
        secret_key.save(sk_stream);
        public_key.save(data_stream);
        relin_keys.save(data_stream);

        // Construction of utility objects.
        // Encryptor encryptor_pk(context, public_key);
        Encryptor encryptor_sk(context, secret_key);
        Evaluator evaluator(context);
        Decryptor decryptor(context, secret_key);
        BatchEncoder batch_encoder(context);
        size_t slot_count = batch_encoder.slot_count();
        uint64_t plain_modulus = get_T_from_context(context);
        uint64_t max_value = sqrt(plain_modulus/2); // ensure the distance won't overflow

        /*
        Generate and encrypt the rider's position data.
        (x, y) in a grid of [0, max_value-1] x [0, max_value-1]
        */
        const uint64_t rider_x = random_uint64() % max_value;
        const uint64_t rider_y = random_uint64() % max_value;
        sk_stream << rider_x << " " << rider_y;

        vector<uint64_t> rider_vector; // rider vector
        for (size_t i = 0; i < slot_count; i++) {
            rider_vector.push_back(i % 2 == 0 ? rider_x : rider_y);
        }
        
        Plaintext rider_plaintext; // encoded plaintext
        batch_encoder.encode(rider_vector, rider_plaintext);
        print_matrix(rider_vector, slot_count/2);

        Serializable<Ciphertext> rider_ciphertext = encryptor_sk.encrypt_symmetric(rider_plaintext);
        rider_ciphertext.save(data_stream);
    }

    /*
    Step 3: Server broadcasts the rider's pk.
    
    Note - This step is empty here because data_stream is shared with drivers as well.
    In reality, server should load pk, rlk and ciphertext of rider from data_stream,
    store locally, and broadcast only pk to the drivers.
    */
    {
    }

    /*
    Step 4: Drivers encrypt their own coordinates with the rider's pk and send
    to the server.
    */
    {
        // Load the parms from the parms_stream and reset the read pointer.
        EncryptionParameters parms;
        parms.load(parms_stream);
        parms_stream.seekg(0, parms_stream.beg);
        SEALContext context(parms);

        // Load rider's pk from data_stream.
        PublicKey public_key;
        public_key.load(context, data_stream);

        // Construction of utility objects.
        Encryptor encryptor_pk(context, public_key);
        BatchEncoder batch_encoder(context);
        size_t slot_count = batch_encoder.slot_count();
        uint64_t plain_modulus = get_T_from_context(context);
        uint64_t max_value = sqrt(plain_modulus/2); // ensure the distance won't overflow

        /*
        Generate and encrypt the drivers' position data. 
        */
        for (size_t i = 0; i < driver_num; ++i) {
            vector<uint64_t> driver_vector(slot_count, 0ULL);
            driver_vector[i<<1] = random_uint64() % max_value;
            driver_vector[(i<<1)+1] = random_uint64() % max_value;
            driver_vectors.push_back(driver_vector);
            if (i < 3) print_matrix(driver_vector, slot_count/2);
            
            Plaintext driver_plaintext;
            batch_encoder.encode(driver_vector, driver_plaintext);

            Serializable<Ciphertext> driver_ciphertext = encryptor_pk.encrypt(driver_plaintext);
            driver_ciphertext.save(data_stream);
        }
    }

    /*
    Step 5: The server evaluates the squared distance between the rider and each of the drivers
    homomorphically with relinearization after each operation and then send back to the rider;
    */
    {
        // Load the parms from the parms_stream and reset the read pointer.
        EncryptionParameters parms;
        parms.load(parms_stream);
        parms_stream.seekg(0, parms_stream.beg);
        SEALContext context(parms);

        // Load rider's rlk from data_stream.
        RelinKeys relin_keys;
        relin_keys.load(context, data_stream);

        // Construction of utility objects.
        Evaluator evaluator(context);

        // Load rider's ciphertext.
        Ciphertext rider_ciphertext;
        rider_ciphertext.load(context, data_stream);

        // Evaluate addition while loading drivers' ciphertexts.
        evaluator.negate_inplace(rider_ciphertext);
        for (size_t i = 0; i < driver_num; ++i) {
            Ciphertext driver_ciphertext;
            driver_ciphertext.load(context, data_stream);
            evaluator.add_inplace(rider_ciphertext, driver_ciphertext);
            evaluator.relinearize_inplace(rider_ciphertext, relin_keys);
        }
        evaluator.square_inplace(rider_ciphertext);
        evaluator.relinearize_inplace(rider_ciphertext, relin_keys);
        rider_ciphertext.save(data_stream);
    }

    /*
    Step 6: The rider decrypts the evaluation result with her sk, and identifies
    the closest driver.
    */
    {   
        // Load the parms from the parms_stream and reset the read pointer.
        EncryptionParameters parms;
        parms.load(parms_stream);
        // parms_stream.seekg(0, parms_stream.beg);
        SEALContext context(parms);

        // Load its own sk and coordinates from sk_stream.
        SecretKey secret_key;
        secret_key.load(context, sk_stream);
        uint64_t rider_x, rider_y;
        sk_stream >> rider_x >> rider_y;

        // Construction of utility objects.
        Evaluator evaluator(context);
        Decryptor decryptor(context, secret_key);
        BatchEncoder batch_encoder(context);
        size_t slot_count = batch_encoder.slot_count();
        uint64_t plain_modulus = get_T_from_context(context);
        uint64_t max_value = sqrt(plain_modulus/2); // ensure the distance won't overflow

        // Load evaluated result ciphertext from data_stream.
        Ciphertext rider_ciphertext;
        rider_ciphertext.load(context, data_stream);
        cout << "    + Noise budget in result: " << decryptor.invariant_noise_budget(rider_ciphertext) << " bits" << endl << endl;

        // Decrypt and decode the result.
        Plaintext result_plaintext;
        vector<uint64_t> result_vector;
        print_line(__LINE__);
        cout << "Decrypt and decode result." << endl;
        decryptor.decrypt(rider_ciphertext, result_plaintext);
        batch_encoder.decode(result_plaintext, result_vector);
        cout << "    + Result plaintext matrix:" << endl;
        print_matrix(result_vector, slot_count/2);

        // Find the closest driver and count the total errors found.
        uint64_t min_id = 0, min_x = max_value, min_y = max_value, min_dist = plain_modulus;
        size_t error_count = 0;

        for (size_t i = 0; i < driver_num; ++i) {
            uint64_t driver_x = driver_vectors[i][i<<1], driver_y = driver_vectors[i][1+(i<<1)];

            uint64_t computed_dist = result_vector[i<<1] + result_vector[1+(i<<1)];
            uint64_t actual_dist = distance(driver_x, driver_y, rider_x, rider_y);
            
            if (computed_dist == actual_dist) {
                if (computed_dist < min_dist) {
                    min_id = i, min_x = driver_x, min_y = driver_y, min_dist = computed_dist;
                }
            }
            else {
                error_count++;
                cout << computed_dist << " " << actual_dist << endl;
            }
        }

        cout << "Rider computes the closest driver position." << endl;
        cout << "Closest driver: driver #" << min_id << " at (" << min_x << ", " << min_y <<") with distance " << min_dist << "." << endl;
        cout << error_count << " errors found." << endl;
    }

    return 0;
}
