function M1_spar = left_mat(miu, m, n)

    count = 1;
    temp = 5 * m * n - 12 * m - 12 * n + 28; % NO. of non-zero entries in the sparse matrix
    temp = zeros(temp, 3);

    for i = 2:m - 3

        for j = 1:n - 2
            index = (i - 1) * (n - 2) + j; % queue of mesh points

            if (j ~= 1) && (j ~= n - 2)% interior

                temp(count, 1) = index;
                temp(count, 2) = index;
                temp(count, 3) = 1 + 2 * miu;
                count = count + 1;

                temp(count, 1) = index;
                temp(count, 2) = index - 1;
                temp(count, 3)=-0.5 * miu;
                count = count + 1;

                temp(count, 1) = index;
                temp(count, 2) = index + 1;
                temp(count, 3)=-0.5 * miu;
                count = count + 1;

                temp(count, 1) = index;
                temp(count, 2) = index - (n - 2);
                temp(count, 3)=-0.5 * miu;
                count = count + 1;

                temp(count, 1) = index;
                temp(count, 2) = index + (n - 2);
                temp(count, 3)=-0.5 * miu;
                count = count + 1;

            elseif (j == 1)

                temp(count, 1) = index;
                temp(count, 2) = index;
                temp(count, 3) = 1 + 1.5 * miu;
                count = count + 1;

                temp(count, 1) = index;
                temp(count, 2) = index + 1;
                temp(count, 3)=-0.5 * miu;
                count = count + 1;

                temp(count, 1) = index;
                temp(count, 2) = index + n - 2;
                temp(count, 3)=-0.5 * miu;
                count = count + 1;

                temp(count, 1) = index;
                temp(count, 2) = index - n + 2;
                temp(count, 3)=-0.5 * miu;
                count = count + 1;

            elseif (j == n - 2)

                temp(count, 1) = index;
                temp(count, 2) = index;
                temp(count, 3) = 1 + 1.5 * miu;
                count = count + 1;

                temp(count, 1) = index;
                temp(count, 2) = index - 1;
                temp(count, 3)=-0.5 * miu;
                count = count + 1;

                temp(count, 1) = index;
                temp(count, 2) = index + n - 2;
                temp(count, 3)=-0.5 * miu;
                count = count + 1;

                temp(count, 1) = index;
                temp(count, 2) = index - n + 2;
                temp(count, 3)=-0.5 * miu;
                count = count + 1;

            end

        end

    end

    temp(count, 1) = 1;
    temp(count, 2) = 1;
    temp(count, 3) = 1 + miu;
    count = count + 1;

    temp(count, 1) = 1;
    temp(count, 2) = 2;
    temp(count, 3)=-0.5 * miu;
    count = count + 1;

    temp(count, 1) = 1;
    temp(count, 2) = n - 1;
    temp(count, 3)=-0.5 * miu;
    count = count + 1;

    temp(count, 1) = n - 2;
    temp(count, 2) = n - 2;
    temp(count, 3) = 1 + miu;
    count = count + 1;

    temp(count, 1) = n - 2;
    temp(count, 2) = n - 3;
    temp(count, 3)=-0.5 * miu;
    count = count + 1;

    temp(count, 1) = n - 2;
    temp(count, 2) = 2 * n - 4;
    temp(count, 3)=-0.5 * miu;
    count = count + 1;

    for j = 2:n - 3
        temp(count, 1) = j;
        temp(count, 2) = j;
        temp(count, 3) = 1 + 1.5 * miu;
        count = count + 1;

        temp(count, 1) = j;
        temp(count, 2) = j - 1;
        temp(count, 3)=-0.5 * miu;
        count = count + 1;

        temp(count, 1) = j;
        temp(count, 2) = j + 1;
        temp(count, 3)=-0.5 * miu;
        count = count + 1;

        temp(count, 1) = j;
        temp(count, 2) = j + n - 2;
        temp(count, 3)=-0.5 * miu;
        count = count + 1;
    end

    temp(count, 1) = (n - 2) * (m - 3) + 1;
    temp(count, 2) = (n - 2) * (m - 3) + 1;
    temp(count, 3) = 1 + miu;
    count = count + 1;

    temp(count, 1) = (n - 2) * (m - 3) + 1;
    temp(count, 2) = (n - 2) * (m - 3) + 2;
    temp(count, 3)=-0.5 * miu;
    count = count + 1;

    temp(count, 1) = (n - 2) * (m - 3) + 1;
    temp(count, 2) = (n - 2) * (m - 3) + 1 - (n - 2);
    temp(count, 3)=-0.5 * miu;
    count = count + 1;

    temp(count, 1) = (n - 2) * (m - 2);
    temp(count, 2) = (n - 2) * (m - 2);
    temp(count, 3) = 1 + miu;
    count = count + 1;

    temp(count, 1) = (n - 2) * (m - 2);
    temp(count, 2) = (n - 2) * (m - 2) - 1;
    temp(count, 3)=-0.5 * miu;
    count = count + 1;

    temp(count, 1) = (n - 2) * (m - 2);
    temp(count, 2) = (n - 2) * (m - 3);
    temp(count, 3)=-0.5 * miu;
    count = count + 1;

    for j = ((n - 2) * (m - 3) + 2):((n - 2) * (m - 2) - 1)
        temp(count, 1) = j;
        temp(count, 2) = j;
        temp(count, 3) = 1 + 1.5 * miu;
        count = count + 1;

        temp(count, 1) = j;
        temp(count, 2) = j - 1;
        temp(count, 3)=-0.5 * miu;
        count = count + 1;

        temp(count, 1) = j;
        temp(count, 2) = j + 1;
        temp(count, 3)=-0.5 * miu;
        count = count + 1;

        temp(count, 1) = j;
        temp(count, 2) = j - n + 2;
        temp(count, 3)=-0.5 * miu;
        count = count + 1;
    end

    temp = temp(1:count - 1, :);

    M1_spar = spconvert(temp);
