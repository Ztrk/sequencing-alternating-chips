#include <algorithm>
#include <iomanip>
#include <iostream>
#include "exact_solver.h"
using namespace std;

#define NOTAVAILABLE 0
#define AVAILABLE 1
#define STOP 0
#define KEEPGOING 1
#define EVEN 0
#define ODD 1

ExactSolver::ExactSolver(const std::vector<std::string> &even_spectrum_input, 
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

    add_oligos(even_spectrum, odd_spectrum);

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
}

void ExactSolver::add_oligos(std::vector<std::string> &even_spectrum, const std::vector<std::string> &odd_spectrum) {
    std::unordered_set<std::string> parts;
    size_t part_len = even_spectrum[0].size() - 2;
    for (std::string &oligo : even_spectrum) {
        for (size_t i = 0; i <= oligo.size() - part_len; i += 2) {
            parts.insert(oligo.substr(i, part_len));
        }
    }

    for (const std::string &oligo : odd_spectrum) {
        std::string part = oligo.substr(0, part_len);
        if (parts.find(part) == parts.end()) {
            even_spectrum.push_back(part);
            parts.insert(part);
        }
    }
}

std::string ExactSolver::solve()
{
    //create even and odd paths starting with given elements
    DNAPath evenPath(evenStart, 0);
    DNAPath oddPath(oddStart, 1);

    //create vertices availability vector
    std::vector <bool> verticesAvailability(even_spectrum.size(), AVAILABLE);

    //set starting vertices as not available
    verticesAvailability[0] = NOTAVAILABLE;
    verticesAvailability[1] = NOTAVAILABLE;

    auto odd_oligos = odd_spectrum;
    int odd_errors = count_errors_odd_spectrum(evenPath, oddPath, oddLength - 1, odd_oligos);

    solveRecursion(oddPath, evenPath, verticesAvailability, {0, odd_errors}, odd_oligos);

    return solutions.size() > 0 ? solutions[0] : "NOT FOUND";
}

bool ExactSolver::solveRecursion(DNAPath shorterPath, DNAPath longerPath, 
        std::vector <bool> verticesAvailability, std::vector <int> errorsCount, unordered_set<string> odd_oligos)
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
    // if all vertices are used
    else if (find(verticesAvailability.begin(), verticesAvailability.end(), AVAILABLE) == verticesAvailability.end())
    {
        //if we were expending longer path swap paths
        //to make longerPath variable actually longer
        //and shorterPath shorter
        if (longerPath.getLength() < shorterPath.getLength()) {
            swap(longerPath, shorterPath);
        }

        return solveRecursionOdd(shorterPath, longerPath, odd_oligos);
    }

    //find all next possible vertices
    std::vector <int> candidates = findAllPossibleVertices(&shorterPath, &longerPath, verticesAvailability);

    //for each vertexID in candidates
    for (int vertexID : candidates)
    {

        std::string newElement = even_spectrum[vertexID];
        
        // Remember previous values
        auto errorsCountPrev = errorsCount;
        DNAPath newPath1 = shorterPath;
        DNAPath newPath2 = longerPath;
        auto odd_oligos_copy = odd_oligos;

        // Add vertex and update even errors
        int even_errors = newPath1.addElement(newElement, vertexID);
        errorsCount[EVEN] += even_errors;

        // When too many even spectrum errors, take next vertex
        if (evenNegativeErrorsExpectedCount < errorsCount[EVEN]) {
            errorsCount[EVEN] -= even_errors;
            continue;
        }

        int odd_errors = count_errors_odd_spectrum(newPath1, newPath2, shorterPath.getLength() + 1, odd_oligos_copy);
        errorsCount[ODD] += odd_errors;

        // Take new oligo if too many errors in odd spectrum
        if (oddNegativeErrorsExpectedCount < errorsCount[ODD]) {
            errorsCount[EVEN] -= even_errors;
            errorsCount[ODD] -= odd_errors;
            continue;
        }

        //set vertex as not available
        verticesAvailability[vertexID] = NOTAVAILABLE;

        //result of new recursion
        bool result;

        //if newPath1 is shorter
        if (newPath1.getLength() <= newPath2.getLength())
        {
            //keep looking for the solutions
            result = solveRecursion(newPath1, newPath2, verticesAvailability, errorsCount, odd_oligos_copy);
        }
        //if newPath2 is shorter
        else if (newPath2.getLength() < newPath1.getLength())
        {
            //keep looking for the solutions
            result = solveRecursion(newPath2, newPath1, verticesAvailability, errorsCount, odd_oligos_copy);
        }

        //set vertex as available
        verticesAvailability[vertexID] = AVAILABLE;
        errorsCount = errorsCountPrev;

        if (result == STOP)
        {
            return STOP;
        }
    }

    // If we are in the ending of sequence
    //do it only when longerPath is actually longer
    //it will prevent from the never ending recursion loop
    if (shorterPath.getLength() + evenLength - 1 > length &&
        shorterPath.getLength() < longerPath.getLength())
    {
        return solveRecursion(longerPath, shorterPath, verticesAvailability, errorsCount, odd_oligos);
    }

    return KEEPGOING;
}

bool ExactSolver::solveRecursionOdd(DNAPath shorterPath, DNAPath longerPath, 
        unordered_set<string> availableOddElements)
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
    //if all odd elements are used
    else if (availableOddElements.size() == 0)
    {
        std::string newSolution = shorterPath.merge(longerPath);
        
        //add solution to solutions vector
        solutions.push_back(newSolution);

        //if found solutions count is equal to solutions limit stop algorithm
        if (static_cast<int>(solutions.size()) == solutionsLimit)
        {
            return STOP;
        }
        
        return KEEPGOING;
    }

    //for each oddElement in availableOddElements
    for (const string &oddElement : availableOddElements)
    {
        //create new paths
        DNAPath newPath1 = shorterPath;
        DNAPath newPath2 = longerPath;

        //add oddElement to new shorter path
        //if adding oddElement was not possible
        //continue to try with next oddElement
        if(!newPath1.addOddElement(longerPath, oddElement))
        {
            continue;
        }

        //delete oddElement from newElements
        auto newElements = availableOddElements;
        count_errors_odd_spectrum(newPath1, newPath2, shorterPath.getLength() + 1, newElements);

        //result of new recursion
        bool result;

        //if newPath1 is shorter
        if (newPath1.getLength() <= newPath2.getLength())
        {
            //keep looking for the solutions
            result = solveRecursionOdd(newPath1, newPath2, newElements);
        }
        //if newPath2 is shorter
        else if (newPath2.getLength() < newPath1.getLength())
        {
            //keep looking for the solutions
            result = solveRecursionOdd(newPath2, newPath1, newElements);
        }

        if (result == STOP)
        {
            return STOP;
        }
    }

    return KEEPGOING;
}

std::vector <int> ExactSolver::findAllPossibleVertices(DNAPath* shorterPath, DNAPath* longerPath, 
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

            return isSecondVerified < isFirstVerified;
        }

        return false;
    };

    //sort outgoingVertices by overlap and isVerified
    std::sort(candidates.begin(), candidates.end(), overlapThenVerificationLambda);

    //return Candidates created with
    //verifiedCandidates and unverifiedCandidates
    return candidates;
}

bool ExactSolver::verifyElementWithOddSpectrum(DNAPath* shorterPath, DNAPath* longerPath, 
        int elementID)
{
    //calculate from which position of dna path the interesting part starts
    //example:
    //shorter path: AXCXAXCXT
    //new element:        CXTXGXT
    //longer path: GXCXGXTXAXGXCX
    //interesting part:  TXAXGG 
    //interesting part length is equal to odd elements length
    if (shorterPath->getLength() > longerPath->getLength()) {
        return true;
    }

    int beginningPosition = shorterPath->getLength() - oddLength + 2;

    //create verifying element to look for in odd spectrum
    //beginning of element is the interesting part described above
    std::string beginning = longerPath->substr(beginningPosition, oddLength - 1);

    //find the expanding nucleotide
    char expandingNucleotide = shorterPath->findExpandingNucleotide(even_spectrum[elementID]);

    //connect it into verifying element
    std::string verifyingElement = beginning + expandingNucleotide;

    // if verifying element is in odd spectrum verification succeeds
    return odd_spectrum.find(verifyingElement) != odd_spectrum.end();
}

int ExactSolver::count_errors_odd_spectrum(DNAPath &extended_path, DNAPath &second_path, 
    int start, unordered_set<string> &odd_oligos)
{
    int errors = 0;
    int end = min(second_path.getLength(), extended_path.getLength());

    // Iterate over every ending of odd oligo
    for (int i = start; i <= end; ++i) {
        // Build oligo
        string verifying_oligo;
        if (i >= second_path.getLength() || second_path.substr(i, 1) == "X") {
            verifying_oligo = second_path.substr(i - oddLength + 1, oddLength - 1);
            verifying_oligo += extended_path.substr(i, 1);
        }
        else {
            verifying_oligo = extended_path.substr(i - oddLength + 1, oddLength - 1);
            verifying_oligo += second_path.substr(i, 1);
        }
        int num_removed = odd_oligos.erase(verifying_oligo);
        if (num_removed == 0) {
            ++errors;
        }
    }

    return errors;
}