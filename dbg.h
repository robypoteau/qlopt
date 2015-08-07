#ifndef __dbg_h__
#define __dbg_h__

#include <iostream>
#include <exception>
#include <string.h>

#ifdef NDEBUG
#define debug(M)
#else
#define debug(M) std::cerr << "DEBUG: "<<  __FILE__ << ":" << __LINE__ << ": "<< M << '\n'
#endif

#define log_err(M) std::cerr << "[ERROR] in "<<  __FILE__ << ":" << __LINE__ << ": "<< M << '\n'
#define log_warn(M) std::cerr << "[WARN] in "<<  __FILE__ << ":" << __LINE__ << ": "<< M << '\n'
#define log_info(M) std::cerr << "[INFO] in "<<  __FILE__ << ":" << __LINE__ << ": "<< M << '\n'

#define note(M) std::cout << "[NOTE] in "<<  __FILE__ << ":" << __LINE__ << ": "<< M << '\n'

#define try_catch(A, M) try { A;} catch(std::exception& e) { log_err(e.what());}

#define check(A, M) if(!(A)) { log_err(M); exit(1); }

#define check_debug(A, M) if(!(A)) { debug(M); }

#define sentinel(M) { log_err(M); exit(1); }

//#define clean(A) if(A) free(A)

//#define clean_file(A) if(A) fclose(A)

#endif
