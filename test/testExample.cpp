// Test example from: https://medium.com/@AlexanderObregon/introduction-to-unit-testing-and-test-driven-development-in-c-1a825a00b1b4

#define BOOST_TEST_MODULE FactorialTest
#include <boost/test/included/unit_test.hpp>

#include <DiscreteElasticRod.hpp>

int Factorial(int n) {
    if (n <= 1) return 1;
    else return n * Factorial(n - 1);
}

// Test case for Factorial function
BOOST_AUTO_TEST_CASE(HandlesZeroInput) {
    BOOST_CHECK_EQUAL(Factorial(0), 1);
}

BOOST_AUTO_TEST_CASE(HandlesPositiveInput) {
    BOOST_CHECK_EQUAL(Factorial(1), 1);
    BOOST_CHECK_EQUAL(Factorial(2), 2);
    BOOST_CHECK_EQUAL(Factorial(3), 6);
    BOOST_CHECK_EQUAL(Factorial(4), 24);
    BOOST_CHECK_EQUAL(Factorial(5), 120);
}
