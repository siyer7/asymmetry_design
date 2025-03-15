% Lookup table to determine instruction prompt used in each trial

function text = getInstructionText(category,trialAxisName,antiTask,promptVariant,equivalentVariantID)

% taskStruct.categoryNames = {'Animals','Cars','Faces','Fruits'};
%taskStruct.axisNames = {{'Colorful','Count'},{'New','Colorful'},...
%    {'New','Identical'},{'Count','Identical'}};

% Singular category names
singularCategoryLabels = {'animal','car','face','item'};
pluralCategoryLabels = {'animals','cars','faces','items'};
alternativeSingularCategoryLabels = {'creature','vehicle','person','object'};
alternativePluralCategoryLabels = {'creatures','vehicles','people','objects'};

if strcmp(trialAxisName,'New')
    if antiTask
        if promptVariant
            if category == 3
                if equivalentVariantID
                    text = ['Avoid the younger ' singularCategoryLabels{category}];
                else
                    text = ['Do not select the youthful ' alternativeSingularCategoryLabels{category}];
                end
            else
                if equivalentVariantID
                    text = ['Avoid the newer ' singularCategoryLabels{category}];
                else
                    text = ['Do not select the more modern ' alternativeSingularCategoryLabels{category}];
                end
            end
        else      
            if category == 3
                if equivalentVariantID
                    text = ['Avoid the older ' singularCategoryLabels{category}];                
                else
                    text = ['Do not select the more aged ' alternativeSingularCategoryLabels{category}];                
                end
            else
                if equivalentVariantID
                    text = ['Avoid the older ' singularCategoryLabels{category}];                
                else
                    text = ['Do not select the less modern ' alternativeSingularCategoryLabels{category}];                
                end
            end
        end
    else
        if promptVariant
            if category == 3
                if equivalentVariantID
                    text = ['Choose the younger ' singularCategoryLabels{category}];
                else
                    text = ['Select the more youthful ' alternativeSingularCategoryLabels{category}];
                end
            else
                if equivalentVariantID
                    text = ['Choose the newer ' singularCategoryLabels{category}];
                else
                    text = ['Select the more modern ' alternativeSingularCategoryLabels{category}];
                end
            end
        else
            if category == 3
                if equivalentVariantID
                    text = ['Choose the older ' singularCategoryLabels{category}];
                else
                    text = ['Select the more aged ' alternativeSingularCategoryLabels{category}];
                end
            else
                if equivalentVariantID
                    text = ['Choose the older ' singularCategoryLabels{category}];
                else
                    text = ['Select the less modern ' alternativeSingularCategoryLabels{category}];
                end
            end
            
        end
    end
elseif strcmp(trialAxisName,'Identical')
    if antiTask
        if promptVariant
            if equivalentVariantID
                text = ['Are these ' pluralCategoryLabels{category} ' different?'];
            else
                text = ['Do these ' alternativePluralCategoryLabels{category} ' have different identities?'];
            end
        else
            if equivalentVariantID
                text = ['Are these different ' pluralCategoryLabels{category} '?'];
            else
                text = ['Are these photos of different ' alternativePluralCategoryLabels{category} '?'];
            end
        end
    else
        if promptVariant
            if equivalentVariantID
                text = ['Are these ' pluralCategoryLabels{category} ' the same?'];
            else
                text = ['Are these photos of the same ' alternativeSingularCategoryLabels{category} '?'];
            end
        else
            if equivalentVariantID
                text = ['Are these ' pluralCategoryLabels{category} ' identical?'];
            else
                text = ['Do these ' alternativePluralCategoryLabels{category} ' have the same identity?'];
            end
        end
    end
elseif strcmp(trialAxisName,'Count')
    if antiTask
        if promptVariant
            if equivalentVariantID
                text = ['Avoid the image with more ' pluralCategoryLabels{category}];
            else
                text = ['Do not select the image with more ' alternativePluralCategoryLabels{category}];
            end
        else
            if equivalentVariantID
                text = ['Avoid the image with fewer ' pluralCategoryLabels{category}];
            else
                text = ['Do not select the image with fewer ' alternativePluralCategoryLabels{category}];
            end
        end
    else
        if promptVariant
            if equivalentVariantID
                text = ['Choose the image with more ' pluralCategoryLabels{category}];
            else
                text = ['Select the image with more ' alternativePluralCategoryLabels{category}];
            end
        else
            if equivalentVariantID
                text = ['Choose the image with fewer ' pluralCategoryLabels{category}];
            else
                text = ['Select the image with fewer ' alternativePluralCategoryLabels{category}];
            end
        end
    end
elseif strcmp(trialAxisName,'Colorful')
    if antiTask
        if promptVariant
            if equivalentVariantID
                text = ['Avoid the more colorful ' singularCategoryLabels{category}];
            else
                text = ['Do not select the ' alternativeSingularCategoryLabels{category} ' with more colors' ];
            end
        else
            if equivalentVariantID
                text = ['Avoid the less colorful ' singularCategoryLabels{category}];
            else
                text = ['Do not select the ' alternativeSingularCategoryLabels{category} ' with fewer colors'];
            end
            
        end
    else
        if promptVariant
            if equivalentVariantID
                text = ['Choose the more colorful ' singularCategoryLabels{category}];
            else
                text = ['Select the ' alternativeSingularCategoryLabels{category} ' with more colors'];
            end
        else
            if equivalentVariantID
                text = ['Choose the less colorful ' singularCategoryLabels{category}];
            else
                text = ['Select the ' alternativeSingularCategoryLabels{category} ' with fewer colors'];
            end            
        end
    end
end
