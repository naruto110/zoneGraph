#!/bin/bash
./inplace 1 1 > ./result/inplace_lj_1 ;
./inplace 1 16 > ./result/inplace_lj_16 ;
./inplace 1 32 > ./result/inplace_lj_32 ;
./inplace 1 64 > ./result/inplace_lj_64 ;

./inplace 9 1 > ./result/inplace_s21_1 ;
./inplace 9 16 > ./result/inplace_s21_16 ;
./inplace 9 32 > ./result/inplace_s21_32 ;
./inplace 9 64 > ./result/inplace_s21_64 ;


