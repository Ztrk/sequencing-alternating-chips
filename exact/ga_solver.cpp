#include <algorithm>
#include <iomanip>
#include <iostream>
#include "ga_solver.h"

#define NOTAVAILABLE 0
#define AVAILABLE 1
#define STOP 0
#define KEEPGOING 1
#define EVEN 0
#define ODD 1

GaSolver::GaSolver(const std::vector<std::string> &even_spectrum_input, 
    const std::vector<std::string> &odd_spectrum, 
    const std::string &start, int length)
: even_spectrum(even_spectrum_input), 
    odd_spectrum(odd_spectrum.begin(), odd_spectrum.end()), length(length) 
{
    std::string even_start = start;
    for (size_t i = 1; i < even_start.size(); i += 2) {
        even_start[i] = 'X';
    }
    std::string odd_start = start.substr(1, start.size() - 2);
    for (size_t i = 1; i < odd_start.size(); i += 2) {
        odd_start[i] = 'X';
    }
    auto start_pos = find(even_spectrum.begin(), even_spectrum.end(), even_start);
    if (start_pos != even_spectrum.end()) {
        even_spectrum.erase(start_pos);
    }
    even_spectrum.insert(even_spectrum.begin(), even_start);
    even_spectrum.insert(even_spectrum.begin() + 1, odd_start);

    //calculate k value
    k = (even_spectrum[0].size() + 1) / 2;

    //calculate evenLength and oddLength values
    evenLength = even_spectrum[0].size();
    oddLength = evenLength - 1;

    //calculate how many even and odd elements 
    //should perfect spectrum have
    perfectSpectrumEvenCount = length - (2 * k) + 2;
    perfectSpectrumOddCount = length - (2 * k) + 3;

    //calculate even and odd negative errors count
    //even_spectrum.size() - 1 because odd path starting element doesn't count
    evenNegativeErrorsExpectedCount = perfectSpectrumEvenCount - (even_spectrum.size() - 1);
    oddNegativeErrorsExpectedCount = perfectSpectrumOddCount - odd_spectrum.size();

    //create overlap between elements graph
    overlapGraph = new OverlapGraph(even_spectrum);

    //find even and odd path starting elements
    //we created them before and put them
    //in the beggining of even_spectrum vector
    evenStart = even_spectrum[0];
    oddStart = even_spectrum[1];

    //save start time
    startTime = std::chrono::high_resolution_clock::now();

    //how many iterations have been rejected
    rejectedItereationsCount = 0;
}

std::string GaSolver::solve()
{
    //create even and odd paths starting with given elements
    DNAPath evenPath(evenStart, 0);
    DNAPath oddPath(oddStart, 1);

    //create vertices availability vector
    std::vector <bool> verticesAvailability(even_spectrum.size(), AVAILABLE);

    //set starting vertices as not available
    verticesAvailability[0] = NOTAVAILABLE;
    verticesAvailability[1] = NOTAVAILABLE;

    solveRecursion(oddPath, evenPath, verticesAvailability, {0, 0});

    return "abc";
}

bool GaSolver::solveRecursion(DNAPath shorterPath, DNAPath longerPath, 
        std::vector <bool> verticesAvailability, std::vector <int> errorsCount)
{
    //calculate alhorithm work duration
    auto duration = std::chrono::high_resolution_clock::now() - startTime;
    //if execution time exceeded the limit stop algorithm
    if (timeLimit < std::chrono::duration_cast <std::chrono::seconds> (duration).count())
    {
        return STOP;
    }
    //if one path is longer than expected length
    else if (length < shorterPath.getLength() ||
        length < longerPath.getLength())
    {
        return KEEPGOING;
    }
    //if dna paths contain enough elements
    else if (longerPath.getLength() == length)
    {
        std::cout << "solution:" << std::endl;
        shorterPath.print();
        longerPath.print();
        std::cout << std::endl;

        //create solution TODO
        std::string newSolution = "abc";

        //test solution TODO
        bool isCorrect = true;
        
        if (isCorrect)
        {
            //add solution to solutions vector
            solutions.push_back(newSolution);

            //if found solutions count is equal to solutions limit stop algorithm
            if (solutions.size() == solutionsLimit)
            {
                return STOP;
            }
        }
        else if (!isCorrect)
        {
            //increment rejected iterations count
            rejectedItereationsCount++;

            //if rejected iterations count exceeded the limit stop algorithm
            if (iterationsLimit < rejectedItereationsCount)
            {
                return STOP;
            }
        }

        return KEEPGOING;
    }

    //find all next possible vertices
    std::vector <int> candidates = findAllPossibleVertices(&shorterPath, &longerPath, verticesAvailability);

    std::cout << "first time. Options:";
    for(int i = 0; i < candidates.size(); i++)
    {
        std::cout << " " << candidates[i]; 
    }
    std::cout << std::endl << std::endl;

    //for each vertexID in candidates
    for (int vertexID : candidates)
    {
        //if vertex is not verified by odd spectrum
        if (!verifyElementWithOddSpectrum(&shorterPath, &longerPath, vertexID))
        {
            //increment odd negative errors count
            errorsCount[ODD]++;
            std::cout << "odd: " << errorsCount[ODD] << std::endl;

            //if odd negative errors count is greater than expected count
            if (oddNegativeErrorsExpectedCount < errorsCount[ODD])
            {
                std::cout << "odd negative" << std::endl << std::endl;
                //decrement odd negative errors count
                errorsCount[ODD]--;

                //take new vertexID
                continue;
            }
        }

        //find vertex element string
        std::string newElement = even_spectrum[vertexID];

        //create new paths
        DNAPath newPath1 = shorterPath;
        DNAPath newPath2 = longerPath;

        //add vertex to new shorter path
        //assumed even negative errors count
        int assumedErrorsCount = newPath1.addElement(newElement, vertexID);
        
        //add assumed errors to error counter
        errorsCount[EVEN] += assumedErrorsCount;
        std::cout << "even: " << errorsCount[EVEN] << std::endl;

        //if even negative errors count is greater than expected count
        if (evenNegativeErrorsExpectedCount < errorsCount[EVEN])
        {
            std::cout << "even negative" << std::endl << std::endl;

            //subtract assumed errors count from even errors counter
            errorsCount[EVEN] -= assumedErrorsCount;

            //take new vertexID
            continue;
        }

        std::cout << "after adding " << vertexID << "  " << newElement << std::endl;
        newPath1.print();
        newPath2.print();
        std::cout << std::endl;

        //set vertex as not available
        verticesAvailability[vertexID] = NOTAVAILABLE;

        //result of new recursion
        bool result;

        //if newPath1 is shorter
        if (newPath1.getLength() <= newPath2.getLength())
        {
            //keep looking for the solutions
            result = solveRecursion(newPath1, newPath2, verticesAvailability, errorsCount);
        }
        //if newPath2 is shorter
        else if (newPath2.getLength() < newPath1.getLength())
        {
            //keep looking for the solutions
            result = solveRecursion(newPath2, newPath1, verticesAvailability, errorsCount);
        }

        //set vertex as available
        verticesAvailability[vertexID] = AVAILABLE;

        std::cout << result << " back to " << newElement << std::endl;
        std::cout << "odd: " << errorsCount[ODD] << std::endl;
        std::cout << "even: " << errorsCount[EVEN] << std::endl;
        shorterPath.print();
        longerPath.print();
        std::cout << std::endl;

        //for (int i = 0; i < verticesAvailability.size(); i++)
        //{
        //    if (verticesAvailability[i] == AVAILABLE)
        //    {
        //        std::cout << i << " ";
        //    }
        //}
        //std::cout << std::endl << std::endl;

        if (result == STOP)
        {
            return STOP;
        }
    }

    return KEEPGOING;
}

std::vector <int> GaSolver::findAllPossibleVertices(DNAPath* shorterPath, DNAPath* longerPath, 
        std::vector <bool> verticesAvailability)
{
    //all outgoing vertices from last vertex
    std::vector <int> outgoingVertices = overlapGraph->getOutgoingVertices(shorterPath->lastElementID);

    //vertices from outgoingVertices that are still available
    std::vector <int> candidates;

    //for each vertex in outgoingVertices
    for (int vertexID : outgoingVertices)
    {
        //if vertex is available add it to candidates
        if (verticesAvailability[vertexID] == AVAILABLE)
        {
            candidates.push_back(vertexID);
        }
    }

    //lambda for getting the better element
    //better means more overlap with the last element
    //if overlaps are equal look for the verification
    //with odd spectrum
    auto overlapThenVerificationLambda = [=](int firstVectorID, int secondVectorID)
    {
        //first element overlap with last element from the path
        int firstOverlap = (overlapGraph->graph)[shorterPath->lastElementID][firstVectorID];

        //second element overlap with last element from the path
        int secondOverlap = (overlapGraph->graph)[shorterPath->lastElementID][secondVectorID];
        
        //elements overlaps differ from each other
        if (firstOverlap != secondOverlap)
        {
            return firstOverlap < secondOverlap;
        }
        //elements overlaps are equal
        else if (firstOverlap == secondOverlap)
        {
            //see if elements are verified by odd spectrum
            //first element verification by odd spectrum
            bool isFirstVerified = verifyElementWithOddSpectrum(shorterPath, longerPath, firstVectorID);

            //second element verification by odd spectrum
            bool isSecondVerified = verifyElementWithOddSpectrum(shorterPath, longerPath, secondVectorID);

            return isSecondVerified <= isFirstVerified;
        }

        return false;
    };

    //sort outgoingVertices by overlap and isVerified
    std::sort(candidates.begin(), candidates.end(), overlapThenVerificationLambda);

    //return Candidates created with
    //verifiedCandidates and unverifiedCandidates
    return candidates;
}

bool GaSolver::verifyElementWithOddSpectrum(DNAPath* shorterPath, DNAPath* longerPath, 
    int elementID)
{
    //calculate from which position of dna path the interesting part starts
    //example:
    //shorter path: AXCXAXCXT
    //new element:        CXTXGXT
    //longer path: GXCXGXTXAXGXCX
    //interesting part:  TXAXGG 
    //interesting part length is equal to odd elements length
    int beginningPosition = shorterPath->getLength() - oddLength + 2;

    //create verifying element to look for in odd spectrum
    //beginning of element is the interesting part described above
    std::string beginning = longerPath->substr(beginningPosition, oddLength - 1);

    //find the expanding nucleotide
    char expandingNucleotide = shorterPath->findExpandingNucleotide(even_spectrum[elementID]);

    //connect it into verifying element
    std::string verifyingElement = beginning + expandingNucleotide;

    //if verifying element is in odd spectrum verification succed
    //else not succed
    if (std::find(odd_spectrum.begin(), odd_spectrum.end(), verifyingElement) != odd_spectrum.end())
    {
        return true;
    }
    else
    {
        return false;
    }
}