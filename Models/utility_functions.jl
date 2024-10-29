# Function to extract theta, alpha, and beta samples and store them in a structured dictionary
function extract_posterior_samples(all_posteriors, sample_vectors, test, start_sample, K)
    n_samples = sum([length(s) for s in sample_vectors])
    n_containers = length(sample_vectors)
    
    current_container = Int[]  # Initialize an empty array to hold the container indices
    current_sample = Int[]     # Initialize an empty array to hold the sample indices

    # Iterate through each container and its samples
    for (container_index, samples) in enumerate(sample_vectors)
        # Append the container index for each sample in the container
        append!(current_container, [container_index for _ in samples])

        # Append the sample index for each sample in the container
        append!(current_sample, [i-1 for i in 1:length(samples)])
    end

    test1 = copy(current_container)
    test2 = copy(current_container)

    #This replaces the last container (test1) with the last container + 1 (test2)
    test2 = append!(test2[1:end-length(sample_vectors[end])], fill(current_container[end]+1, length(sample_vectors[end])))

    if test == "test1"
        current_container = test1
    elseif test == "test2"
        current_container = test2
    end
    

    # Create a nested dictionary to store posterior samples by container and sample index
    posterior_samples = Dict{Int, Dict{Symbol, Vector{Vector{Float64}}}}()

    for i in start_sample:length(current_sample)
        container_index = current_container[i]
        sample_index = current_sample[i]

        # Extract posterior samples for the current container and sample
        theta_samples = vec(all_posteriors[i]["theta[$container_index][$(K)]"])
        alpha_samples = vec(all_posteriors[i]["alpha"])
        beta_samples = vec(all_posteriors[i]["beta[$(K)]"])

        # Initialize nested dictionary if it doesn't exist
        if !haskey(posterior_samples, container_index)
            posterior_samples[container_index] = Dict(
                :theta => Vector{Vector{Float64}}(),
                :alpha => Vector{Vector{Float64}}(),
                :beta => Vector{Vector{Float64}}()
            )
        end

        # Store the samples in the respective vectors
        push!(posterior_samples[container_index][:theta], theta_samples)
        push!(posterior_samples[container_index][:alpha], alpha_samples)
        push!(posterior_samples[container_index][:beta], beta_samples)
    end

    return posterior_samples
end