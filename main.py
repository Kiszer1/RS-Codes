import numpy as np
from sage.all import *

# *************** Utilities ***************

field = 257


# Transforms a message into an ascii array.
def word_to_array(word):
    arr = []
    for ch in word:
        arr.append(ord(ch))
    return arr


# Transforms a polynomial into a message using its coefficients
def poly_to_msg(poly):
    msg = ""
    for coeff in poly.coefficients():
        msg = chr(coeff) + msg
    return msg


# Inserting random numbers up to the field in errors different locations
# of the encoded message
def insert_errors(c, errors):
    # using a set to avoid picking the same location
    locations = set()
    for i in range(errors):
        location = randrange(len(c))
        # if already chose this location
        while location in locations:
            location = randrange(len(c))
        locations.add(location)
        error = randrange(field)
        # making sure our random number is different then the original
        if error == c[location][1]:
            c[location][1] = (c[location] + 1) % field
        else:
            c[location][1] = error
    print(f"Received message = {c}\n")


# *************** RS Encoder ***************

# Creating a polynomial using the ascii value of the message characters as coefficients
def message_to_poly(message):
    poly = 0
    x = PolynomialRing(GF(field), 1, 'x').gen()
    for i in range(len(message)):
        poly += message[i] * x ^ i
    return poly


def rs_encoder(a, message):
    # creating the encoded polynomial
    message_poly = message_to_poly(message)
    # plugging the values of a in the polynomial to receive an encoded message
    return [[element, message_poly(element)] for element in a]


# *************** Unique Decoding ***************

# Building a matrix from the equations
# y_i * E(a_i) = N(a_i) 1 <= i <= n
def build_matrix(a, c, n, k, e):
    matrix = []
    b = []
    for i in range(n):
        matrix.append([])
        for j in range(e):
            matrix[i].append((c[i] * (a[i] ** j)) % field)
        for j in range(e + k):
            matrix[i].append(-(a[i] ** j) % field)
        b.append(-(c[i] * (a[i] ** 2)) % field)
    return b, matrix


# Generating 2 polynomials from a list of coefficients
# E(x) monic of degree e
# N(x) of degree e + k - 1
def make_polynomials(coefficients, e):
    x = PolynomialRing(GF(field), 1, 'x').gen()
    n_poly, e_poly = 0, 0
    for i in range(e):
        e_poly += coefficients[i] * x ** i
    e_poly += x ** e
    for i in range(len(coefficients) - e):
        n_poly += coefficients[i + e] * x ** i
    return e_poly, n_poly


def unique_decoding(n, k, e, encoded_message):
    print(f"\nUnique Decoding:\n")
    # checking for e > 0
    if e <= 0:
        print(f"e = {e}, but e must be greater then 0\n")
        return
    a, y = [], []
    for pair in encoded_message:
        a.append(pair[0])
        y.append(pair[1])
    M = MatrixSpace(GF(field), n, 2 * e + k)
    b, matrix = build_matrix(a, y, n, k, e)
    A = M(matrix)
    try:
        # Solving for x in A * x = b to get E(x) and N(x) coefficients
        coefficients = A.solve_right(vector(b))
    except:
        return
    e_poly, n_poly = make_polynomials(coefficients, e)
    # P(X) <- N(X) / E(X)
    p_x = n_poly // e_poly
    errors = 0
    for i in range(n):
        if p_x(a[i]) != y[i]:
            errors += 1
    # check number of errors in the received p_x
    if errors > e:
        print(f"Failed, got Polynomial = {p_x}\n")
        return
    print(f"Message Polynomial : {p_x}\n")
    print(f"Error corrected message : {poly_to_msg(p_x)}\n")
    return p_x


# *************** List Decoding ***************

# Building a matrix from the equations
# Q(a_i, y_i) = 0,  1 < = i <= n
def build_matrix2(a, c, n, deg_x, deg_y):
    matrix = []
    for i in range(n):
        matrix.append([])
        for j in range(deg_x + 1):
            for r in range(deg_y + 1):
                matrix[i].append((((a[i] ** j) % field) * ((c[i] ** r) % field)) % field)
    return matrix


# Generating a bi-variate polynomial of with deg_x and deg_y from a list of coefficients
# Using Q(X, Y) = c_ij*X**i * y**j , 0 <= i <= deg_x, 0 <= j <= deg_y
def build_bi_polynomial(coefficients, deg_x, deg_y):
    polynomial = 0
    x, y = PolynomialRing(GF(field), 2, ['x', 'y']).gens()
    for i in range(deg_x + 1):
        for j in range(deg_y + 1):
            polynomial += coefficients[((deg_y + 1) * i) + j] * x ** i * y ** j
    return polynomial


def build_q(matrix, n, deg_x, deg_y):
    M = MatrixSpace(GF(field), n, (deg_x + 1) * (deg_y + 1))
    A = M(matrix)
    # Solving For x in A * x = 0
    # Will generate multiple answers
    Qs = A.right_kernel_matrix(basis="computed")
    Q = []
    # Generating a polynomial for each of the answers
    for coefficients in Qs:
        Q.append(build_bi_polynomial(coefficients, deg_x, deg_y))
    return Q


# Checking the factor is in the form Y - P(x)
# Additionally will check that the number of errors in P(x) < e
def good_factor(factor, k, c, e):
    for degrees in factor.dict():
        if degrees[1] > 1 or (degrees[0] > 0 and degrees[1] >= 1) or degrees[0] >= k:
            return False
    errors = 0
    y = PolynomialRing(GF(field), 1, 'y').gen()
    new_factor = y - factor
    for i in range(len(c)):
        if new_factor(i, 0) != c[i]:
            errors += 1
    if errors > e:
        return False
    return True


# Testing all factors of Q(X, Y)
def good_factors(qs, k, c, e):
    L = []
    y = PolynomialRing(GF(field), 1, 'y').gen()
    for q in qs:
        for factor in q.factor():
            if good_factor(factor[0], k, c, e):
                L.append(y - factor[0])
    return L


def list_decoding1(n, k, e, y):
    deg_x = math.floor(math.sqrt(n * (k - 1)))
    # In original code l > 0, but found out this works instead and now we can use l = 0
    if deg_x == 0:
        deg_y = 1
    else:
        deg_y = math.floor(n / math.sqrt(n * (k - 1)))
    a = []
    c = []
    for pair in y:
        a.append(pair[0])
        c.append(pair[1])
    matrix = np.array(build_matrix2(a, c, n, deg_x, deg_y))
    qs = build_q(matrix, n, deg_x, deg_y)
    # Don't want duplicate factors
    L = set(good_factors(qs, k, c, e))

    print(f"\nList Decoding :\n")
    print(f"L :\n{L}\n")
    for poly in L:
        print(f"Error corrected message : {poly_to_msg(poly)}\n")
    return L


def menu():
    input_message = input("Message to send : ")
    k = len(input_message)
    print(f"k = {k}\n")
    e = int(input(f"Number of errors : "))
    n = int(input(f"N : "))
    print(f"e = {e}, n = {n}\n")
    a = list(range(0, n))
    print(f"a : {a}\n")
    message = word_to_array(input_message)
    print(f"Message before encoding : {message}\n")
    encoded_message = rs_encoder(a, message)
    print(f"Encoded message = {encoded_message}\n")
    insert_errors(encoded_message, e)

    unique_decoding(n, k, e, encoded_message)
    list_decoding1(n, k, e, encoded_message)


# *************** Tests ***************
def run(message, e, n):
    a = list(range(n))
    k = len(message)
    encoded_message = rs_encoder(a, word_to_array(message))
    insert_errors(encoded_message, e)
    unique_decoding(n, k, e, encoded_message)
    list_decoding1(n, k, e, encoded_message)


def test1():
    print("\nTest 1\n")
    message = "a"
    e = 0
    n = 1
    run(message, e, n)


def test2():
    print("\nTest 2\n")
    message = "a"
    e = 1
    n = 1
    run(message, e, n)


def test3():
    print("\nTest 3\n")
    message = "ab"
    e = 0
    n = 2
    run(message, e, n)


def test4():
    print("\nTest 4\n")
    message = "ab"
    e = 1
    n = 4
    run(message, e, n)


def test5():
    print("\nTest 5\n")
    message = "ab"
    e = 2
    n = 5
    run(message, e, n)


def test6():
    print("\nTest 6\n")
    message = "abc"
    e = 2
    n = 7
    run(message, e, n)


def test7():
    print("\nTest 7\n")
    message = "abc"
    e = 2
    n = 8
    run(message, e, n)


def test8():
    print("\nTest 8\n")
    message = "abc"
    e = 3
    n = 8
    run(message, e, n)


def test9():
    print("\nTest 9\n")
    message = "abc"
    e = 3
    n = 10
    run(message, e, n)


def test10():
    print("\nTest 10\n")
    message = "abc"
    e = 4
    n = 10
    run(message, e, n)


def test11():
    print("\nTest 11\n")
    message = "abc"
    e = 17
    n = 30
    run(message, e, n)


def test12():
    print("\nTest 12\n")
    message = "abc"
    e = 18
    n = 30
    run(message, e, n)


# TODO
def test13():
    print("\nTest 13\n")
    message = "a"
    e = 0
    n = 1
    run(message, e, n)


def test14():
    print("\nTest 14\n")
    message = "a"
    e = 1
    n = 1
    run(message, e, n)


def test15():
    print("\nTest 15\n")
    message = "ab"
    e = 0
    n = 2
    run(message, e, n)


def test16():
    print("\nTest 16\n")
    message = "ab"
    e = 1
    n = 4
    run(message, e, n)


def test17():
    print("\nTest 17\n")
    message = "ab"
    e = 2
    n = 5
    run(message, e, n)


def test18():
    print("\nTest 18\n")
    message = "abc"
    e = 2
    n = 7
    run(message, e, n)


def test19():
    print("\nTest 19\n")
    message = "abc"
    e = 2
    n = 8
    run(message, e, n)


def test20():
    print("\nTest 20\n")
    message = "abc"
    e = 3
    n = 8
    run(message, e, n)


def test21():
    print("\nTest 21\n")
    message = "abc"
    e = 3
    n = 10
    run(message, e, n)


def test22():
    print("\nTest 22\n")
    message = "abc"
    e = 4
    n = 10
    run(message, e, n)


def test23():
    print("\nTest 23\n")
    message = "abc"
    e = 15
    n = 30
    run(message, e, n)


def test24():
    print("\nTest 24\n")
    message = "abc"
    e = 18
    n = 30
    run(message, e, n)


def tests():
    test1()
    test2()
    test3()
    test4()
    test5()
    test6()
    test7()
    test8()
    test9()
    test10()
    test11()
    test12()
    test13()
    test14()
    test15()
    test16()
    test17()
    test18()
    test19()
    test20()
    test21()
    test22()
    test23()
    test24()


if __name__ == '__main__':
    menu()
