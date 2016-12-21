function OUTPUT_FILE_PATHS = run_piv_job_list(JOBLIST_INPUT)

% First step: Verify that all files that the job
% refers to can be located on the current filesystem. 
% % % % WRITE THIS % % % % %%
% Things to test: 
verify_job_list_file_paths(JOBLIST_INPUT);
% 
% 1) Check existences of all files
%   - Flow images
%   - Mask images
%   - Vector fields for iterative methods
% 2) Check compatibility of iterative methods?
% % % % WRITE THIS % % % % %

% Copy the Job List
JOBLIST_OUTPUT = JOBLIST_INPUT;

% Count the number of jobs
num_jobs = length(JOBLIST_INPUT);

% Loop over all the jobs
for n = 1 : num_jobs
    
    % Extract the job file
    JobFile = JOBLIST_INPUT(n);
    
    % Save the original job file 
    % to the structure before any
    % modifications happen.
    JobList_output(n).InputJobFile = JobFile;
        
    % Run the PIV job file
     OUTPUT_FILE_PATHS{n} = run_piv_job_file(JobFile);
    
end

end