# RS-Codes
An implementation of Reed Solomon's encoding and Unique \ List Decoding
The implementation is based on chapters 5, 11 and 12 of the book Essential Coding Theory by Madhu Sudan.

https://cse.buffalo.edu/faculty/atri/courses/coding-theory/book/web-coding-book.pdf

# Introduction
Reed-Solomon (RS) codes are a type of error-correcting code that was first introduced by Irving S. Reed and Gustave Solomon in 1960. They are widely used in various communication and storage applications, such as CDs, DVDs, QR codes, and satellite communication systems.
RS codes are non-binary linear block codes that work with symbols of any size, such as bytes or groups of bytes. They are capable of correcting burst errors, which occur when multiple errors occur in a row. The Reed-Solomon code is constructed by selecting a finite field and choosing a generator polynomial of degree n-k, where n is the codeword length and k is the number of input symbols.
In our implementation we chose a finite field of 257 because we are using ASCII which has 256 values and 257 is the nearest prime number.

# How to run the file
1. use Jupiter notebook.
2. open the file with the Jupiter notebook.
3. run the tests included.
