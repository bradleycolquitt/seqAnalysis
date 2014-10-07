### Define linear Gaussian system 
### Generate random data set

function random_gauss(n)
    data = zeros((n, 3))


    for i in 1:size(data, 1)
        data[i,1] = randn()
        data[i,2] = 2*data[i,1] + randn()
        data[i,3] = 4*data[i,2] + randn()
    end
    return data
end    




#writedlm("gauss_data.csv", data, ",")
