#include "SignalManager.h"
#include <cstdio>
#include <cstdlib>
#ifdef SIGNAL

//#include "ARBDException.h"

void SignalManager::segfault_handler(int sig, siginfo_t *info, void *secret) 
{
    //comment out for now
    //if(sig != SIGSEGV) 
        //throw ARBDException("segfault_handler should handle segmentation faults only, got %d", sig);
    void *trace[16];
    char **messages = (char **) NULL;
    int i, trace_size = 0;
    ucontext_t *uc = (ucontext_t *) secret;
    //write to stdout for now
    fprintf(stdout, "Segmentation fault identified, faulty address is %p, from %p", info->si_addr, (void *) uc->uc_mcontext.gregs[MY_REG_RIP]);

    trace_size = backtrace(trace, 16);
    /* overwrite sigaction with caller's address */
    trace[1] = (void *) uc->uc_mcontext.gregs[MY_REG_RIP];

    messages = backtrace_symbols(trace, trace_size);
    /* skip first stack frame (points here) */
    
    fprintf(stdout, "Execution path:");
    for (i = 1; i < trace_size; ++i) 
        fprintf(stdout, "\t[bt] %s\n", messages[i]);
    //throw here or just exit ?
    exit(0);
}

void SignalManager::manage_segfault() 
{
	struct sigaction sa;

	sa.sa_sigaction = segfault_handler;
	sigemptyset(&sa.sa_mask);
	sa.sa_flags = SA_RESTART | SA_SIGINFO;

	sigaction(SIGSEGV, &sa, NULL);
}

#else
void SignalManager::segfault_handler(int sig, siginfo_t *info, void *secret) {}
void SignalManager::manage_segfault() {}

#endif
