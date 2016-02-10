#!/usr/bin/sbcl --script
;run like : ./runall.lisp <all args you need except -l>, runs one job per runlist
(map nil #'(lambda (x)
	     (sb-ext:run-program "./tbmon" (concatenate 'list (cdr sb-ext:*posix-argv*) 
							(list "-l" (concatenate 'string "runlists/" x)))
				 :wait t :output t))
     (list "Oct09_H8_Bat_FieldON.0degree"
	   "Oct09_H8_Bat_FieldON.minus15degree"
	   "Oct09_H8_Bat_FieldON.minus22degree"
	   "Oct09_H8_Bat_FieldON.minus30degree"
	   "Oct09_H8_Bat_FieldON.minus7degree"
	   "Oct09_H8_Bat_FieldON.plus22degree"
	   "Oct09_H8_Bat_FieldON.plus30degree"
	   "Oct09_H8_Bat_FieldON.plus7degree"))
