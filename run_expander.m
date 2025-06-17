
% save_file = [mfilename,'_linear_syn.mat']; % name of the file to store results, by default the same as the filename of this script
flag_visual = true; 
mcmc_batch_size = 500;

while true % expand forever unless killed manually
    
    %here is a chain size stop criteria
    if chain.i(end) >= 10000
        break
    end

    chain = chainer_main(chain,[],[],true,[]); % resets chain's acceptance record
    
    chain = chainer_main(chain,mcmc_batch_size,[],true,flag_visual); % expands chain's history
    
    if chain.sizeGB > 1.0
        % thin chain's history, if chain is too big
        chain = chainer_main(chain,-fix(chain.length/2),[],true,[],[]); 
    end
    
    % save chain for future processing
    save(save_file,'chain','save_file') 
    disp(['SAVED: ', save_file])


end    
