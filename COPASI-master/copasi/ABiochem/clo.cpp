/* Begin CVS Header
   $Source: /Volumes/Home/Users/shoops/cvs/copasi_dev/copasi/ABiochem/clo.cpp,v $
   $Revision: 1.7 $
   $Name:  $
   $Author: shoops $
   $Date: 2006/04/27 01:33:58 $
   End CVS Header */

// Copyright � 2005 by Pedro Mendes, Virginia Tech Intellectual
// Properties, Inc. and EML Research, gGmbH.
// All rights reserved.

/*
 * This file was automaticly generated by clo++ version 0.6.1 (REL-0_6_1)
 * http://pmade.org/pjones/software/clo++/
 */

/** @file
 * This file contains the implementation of the parser class.
 * You can control the name of this source file using the
 * cxx_source_file variable.
 **/

#include "clo.h"
#include <cstring>
#include <cstdlib>
#include <cctype>

namespace
  {
  const char const_usage[] =
    "  --specific-activation     use specific activation kinetics\n"
    "  --version                 output version information and exit\n"
    "  -I, --interaction double  the probability of non-additive pairs of inputs\n"
    "  -c, --coop double         the value of Hill coefficients\n"
    "  -g, --genes int           the total number of genes\n"
    "  -i, --inputs int          the number of total interactions\n"
    "  -k, --constants double    the value of inhibition/activation constants\n"
    "  -p, --positive double     the probability of positive inputs\n"
    "  -r, --rates double        the value of rate constants\n"
    "  -s, --seed int            a seed for the random number generator\n"
    "  -t, --total int           the total number of networks to generate\n"
    "  -w, --rewire double       the probability of rewiring links at random\n"
    "  -x, --prefix string       Prefix for filenames\n";

  const char const_help_comment[] =
    "use the -h option for help";

  const char* expand_long_name (const std::string &name);
}

//#########################################################################
clo::parser::parser (void)
    : state_(state_option)
{
  memset(&locations_, 0, sizeof(locations_));
}
//#########################################################################
void clo::parser::parse (int argc, char *argv[], bool call_finalize)
{
  for (int i = 1; i < argc; ++i) parse_element(argv[i], i, source_cl);
  if (call_finalize) finalize();
}
//#########################################################################
void clo::parser::finalize (void)
{
  if (state_ == state_value)
    {
      switch (openum_)
        {
        case option_constants:
          throw option_error("missing value for 'constants' option");
        case option_coop:
          throw option_error("missing value for 'coop' option");
        case option_genes:
          throw option_error("missing value for 'genes' option");
        case option_inputs:
          throw option_error("missing value for 'inputs' option");
        case option_interaction:
          throw option_error("missing value for 'interaction' option");
        case option_positive:
          throw option_error("missing value for 'positive' option");
        case option_prefix:
          throw option_error("missing value for 'prefix' option");
        case option_rates:
          throw option_error("missing value for 'rates' option");
        case option_rewire:
          throw option_error("missing value for 'rewire' option");
        case option_seed:
          throw option_error("missing value for 'seed' option");
        case option_specific_activation:
          throw option_error("missing value for 'specific-activation' option");
        case option_total:
          throw option_error("missing value for 'total' option");
        case option_version:
          throw option_error("missing value for 'version' option");
        }
    }
}
//#########################################################################
void clo::parser::parse_element (const char *element, int position, opsource source)
{
  size_t length = strlen(element);

  switch (state_)
    {
    case state_consume:
      non_options_.push_back(element);
      break;
    case state_option:
      if (length >= 2 && element[0] == '-' && element[1] == '-')
        {
          if (length == 2) {state_ = state_consume; return;}
          element += 2;
          const char *value = element;
          while (*value != 0 && *value != '=') ++value;
          if (*value == '=')
            {
              std::string selement(element, value - element), svalue(++value);

              parse_long_option(selement.c_str(), position, source);
              if (state_ != state_value)
                {
                  std::string error("the '"); error += element; error += "' option does not take a value";
                  throw option_error(error);
                }

              parse_value(svalue.c_str());
              state_ = state_option;
            }
          else
            {
              parse_long_option(element, position, source);
            }
        }
      else if (length >= 2 && element[0] == '-')
        {
          ++element;

          if (length > 2)
            {
              while (*element != 0)
                {
                  parse_short_option(*element, position, source);
                  ++element;

                  if (state_ == state_value && *element != 0)
                    {
                      parse_value(element);
                      state_ = state_option;
                      break;
                    }
                }
            }
          else
            {
              parse_short_option(*element, position, source);
            }
        }
      else
        {
          non_options_.push_back(element);
        }
      break;
    case state_value:
      parse_value(element);
      state_ = state_option;
      break;
    }
}
//#########################################################################
void clo::parser::parse_short_option (char option, int position, opsource source)
{
  switch (option)
    {
    case 'I':
      openum_ = option_interaction;
      state_ = state_value;
      locations_.interaction = position;
      return;
    case 'c':
      openum_ = option_coop;
      state_ = state_value;
      locations_.coop = position;
      return;
    case 'g':
      openum_ = option_genes;
      state_ = state_value;
      locations_.genes = position;
      return;
    case 'i':
      openum_ = option_inputs;
      state_ = state_value;
      locations_.inputs = position;
      return;
    case 'k':
      openum_ = option_constants;
      state_ = state_value;
      locations_.constants = position;
      return;
    case 'p':
      openum_ = option_positive;
      state_ = state_value;
      locations_.positive = position;
      return;
    case 'r':
      openum_ = option_rates;
      state_ = state_value;
      locations_.rates = position;
      return;
    case 's':
      openum_ = option_seed;
      state_ = state_value;
      locations_.seed = position;
      return;
    case 't':
      openum_ = option_total;
      state_ = state_value;
      locations_.total = position;
      return;
    case 'w':
      openum_ = option_rewire;
      state_ = state_value;
      locations_.rewire = position;
      return;
    case 'x':
      openum_ = option_prefix;
      state_ = state_value;
      locations_.prefix = position;
      return;
    case 'h':
      if (source != source_cl) break;
      throw autoexcept(autothrow_help, const_usage);
    case '?':
      if (source != source_cl) break;
      throw autoexcept(autothrow_help, const_usage);
    }

  std::string error("unknown option: '"); error += option; error += "'";
  throw option_error(error);
}
//#########################################################################
void clo::parser::parse_long_option (const char *option, int position, opsource source)
{
  option = expand_long_name(option);

  if (strcmp(option, "constants") == 0)
    {
      openum_ = option_constants;
      locations_.constants = position;
      state_ = state_value;
      return;
    }
  else if (strcmp(option, "coop") == 0)
    {
      openum_ = option_coop;
      locations_.coop = position;
      state_ = state_value;
      return;
    }
  else if (strcmp(option, "genes") == 0)
    {
      openum_ = option_genes;
      locations_.genes = position;
      state_ = state_value;
      return;
    }
  else if (strcmp(option, "inputs") == 0)
    {
      openum_ = option_inputs;
      locations_.inputs = position;
      state_ = state_value;
      return;
    }
  else if (strcmp(option, "interaction") == 0)
    {
      openum_ = option_interaction;
      locations_.interaction = position;
      state_ = state_value;
      return;
    }
  else if (strcmp(option, "positive") == 0)
    {
      openum_ = option_positive;
      locations_.positive = position;
      state_ = state_value;
      return;
    }
  else if (strcmp(option, "prefix") == 0)
    {
      openum_ = option_prefix;
      locations_.prefix = position;
      state_ = state_value;
      return;
    }
  else if (strcmp(option, "rates") == 0)
    {
      openum_ = option_rates;
      locations_.rates = position;
      state_ = state_value;
      return;
    }
  else if (strcmp(option, "rewire") == 0)
    {
      openum_ = option_rewire;
      locations_.rewire = position;
      state_ = state_value;
      return;
    }
  else if (strcmp(option, "seed") == 0)
    {
      openum_ = option_seed;
      locations_.seed = position;
      state_ = state_value;
      return;
    }
  else if (strcmp(option, "specific-activation") == 0)
    {
      openum_ = option_specific_activation;
      locations_.specific_activation = position;
      options_.specific_activation = !options_.specific_activation;
      return;
    }
  else if (strcmp(option, "total") == 0)
    {
      openum_ = option_total;
      locations_.total = position;
      state_ = state_value;
      return;
    }
  else if (strcmp(option, "version") == 0)
    {
      openum_ = option_version;
      locations_.version = position;
      options_.version = !options_.version;
      return;
    }
  else if (source == source_cl && strcmp(option, "help") == 0)
    {
      throw autoexcept(autothrow_help, const_usage);
    }

  std::string error("unknown option '"); error += option; error += "'";
  throw option_error(error);
}
//#########################################################################
void clo::parser::parse_value (const char *value)
{
  switch (openum_)
    {
    case option_constants:
      {
        char *endptr; double tmp = strtod(value, &endptr);
        while (*endptr != 0 && isspace(*endptr)) ++endptr;

        if (*endptr != 0)
          {
            std::string error("invalid floating point value '"); error += value; error += "'";
            throw option_error(error);
          }
        if (tmp < 0.0)
          {
            throw option_error("floating point value out of range, 'constants' min is 0.0");
          }
        options_.constants = tmp;
      }
      break;
    case option_coop:
      {
        char *endptr; double tmp = strtod(value, &endptr);
        while (*endptr != 0 && isspace(*endptr)) ++endptr;

        if (*endptr != 0)
          {
            std::string error("invalid floating point value '"); error += value; error += "'";
            throw option_error(error);
          }
        if (tmp < 0.0)
          {
            throw option_error("floating point value out of range, 'coop' min is 0.0");
          }
        options_.coop = tmp;
      }
      break;
    case option_genes:
      {
        char *endptr; int tmp = strtol(value, &endptr, 0);
        while (*endptr != 0 && isspace(*endptr)) ++endptr;

        if (*endptr != 0)
          {
            std::string error("invalid integer value '"); error += value; error += "'";
            throw option_error(error);
          }
        if (tmp < 1)
          {
            throw option_error("integer value out of range, 'genes' min is 1");
          }
        options_.genes = tmp;
      }
      break;
    case option_inputs:
      {
        char *endptr; int tmp = strtol(value, &endptr, 0);
        while (*endptr != 0 && isspace(*endptr)) ++endptr;

        if (*endptr != 0)
          {
            std::string error("invalid integer value '"); error += value; error += "'";
            throw option_error(error);
          }
        if (tmp < 1)
          {
            throw option_error("integer value out of range, 'inputs' min is 1");
          }
        options_.inputs = tmp;
      }
      break;
    case option_interaction:
      {
        char *endptr; double tmp = strtod(value, &endptr);
        while (*endptr != 0 && isspace(*endptr)) ++endptr;

        if (*endptr != 0)
          {
            std::string error("invalid floating point value '"); error += value; error += "'";
            throw option_error(error);
          }
        if (tmp < 0.0)
          {
            throw option_error("floating point value out of range, 'interaction' min is 0.0");
          }
        if (tmp > 1.0)
          {
            throw option_error("floating point value out of range, 'interaction' max is 1.0");
          }
        options_.interaction = tmp;
      }
      break;
    case option_positive:
      {
        char *endptr; double tmp = strtod(value, &endptr);
        while (*endptr != 0 && isspace(*endptr)) ++endptr;

        if (*endptr != 0)
          {
            std::string error("invalid floating point value '"); error += value; error += "'";
            throw option_error(error);
          }
        if (tmp < 0.0)
          {
            throw option_error("floating point value out of range, 'positive' min is 0.0");
          }
        if (tmp > 1.0)
          {
            throw option_error("floating point value out of range, 'positive' max is 1.0");
          }
        options_.positive = tmp;
      }
      break;
    case option_prefix:
      {
        options_.prefix = value;
      }
      break;
    case option_rates:
      {
        char *endptr; double tmp = strtod(value, &endptr);
        while (*endptr != 0 && isspace(*endptr)) ++endptr;

        if (*endptr != 0)
          {
            std::string error("invalid floating point value '"); error += value; error += "'";
            throw option_error(error);
          }
        if (tmp < 0.0)
          {
            throw option_error("floating point value out of range, 'rates' min is 0.0");
          }
        options_.rates = tmp;
      }
      break;
    case option_rewire:
      {
        char *endptr; double tmp = strtod(value, &endptr);
        while (*endptr != 0 && isspace(*endptr)) ++endptr;

        if (*endptr != 0)
          {
            std::string error("invalid floating point value '"); error += value; error += "'";
            throw option_error(error);
          }
        if (tmp < 0.0)
          {
            throw option_error("floating point value out of range, 'rewire' min is 0.0");
          }
        if (tmp > 1.0)
          {
            throw option_error("floating point value out of range, 'rewire' max is 1.0");
          }
        options_.rewire = tmp;
      }
      break;
    case option_seed:
      {
        char *endptr; int tmp = strtol(value, &endptr, 0);
        while (*endptr != 0 && isspace(*endptr)) ++endptr;

        if (*endptr != 0)
          {
            std::string error("invalid integer value '"); error += value; error += "'";
            throw option_error(error);
          }
        options_.seed = tmp;
      }
      break;
    case option_specific_activation:
      break;
    case option_total:
      {
        char *endptr; int tmp = strtol(value, &endptr, 0);
        while (*endptr != 0 && isspace(*endptr)) ++endptr;

        if (*endptr != 0)
          {
            std::string error("invalid integer value '"); error += value; error += "'";
            throw option_error(error);
          }
        if (tmp < 1)
          {
            throw option_error("integer value out of range, 'total' min is 1");
          }
        options_.total = tmp;
      }
      break;
    case option_version:
      break;
    }
}
//#########################################################################
const char* clo::option_error::get_help_comment (void) const
  {
    return const_help_comment;
  }
//#########################################################################
namespace
  {
  const char* expand_long_name (const std::string &name)
  {
    std::string::size_type name_size = name.size();
    std::vector<const char*> matches;

    if (name_size <= 9 && name.compare(0, name_size, "constants", name_size) == 0)
      matches.push_back("constants");

    if (name_size <= 4 && name.compare(0, name_size, "coop", name_size) == 0)
      matches.push_back("coop");

    if (name_size <= 5 && name.compare(0, name_size, "genes", name_size) == 0)
      matches.push_back("genes");

    if (name_size <= 6 && name.compare(0, name_size, "inputs", name_size) == 0)
      matches.push_back("inputs");

    if (name_size <= 11 && name.compare(0, name_size, "interaction", name_size) == 0)
      matches.push_back("interaction");

    if (name_size <= 8 && name.compare(0, name_size, "positive", name_size) == 0)
      matches.push_back("positive");

    if (name_size <= 6 && name.compare(0, name_size, "prefix", name_size) == 0)
      matches.push_back("prefix");

    if (name_size <= 5 && name.compare(0, name_size, "rates", name_size) == 0)
      matches.push_back("rates");

    if (name_size <= 6 && name.compare(0, name_size, "rewire", name_size) == 0)
      matches.push_back("rewire");

    if (name_size <= 4 && name.compare(0, name_size, "seed", name_size) == 0)
      matches.push_back("seed");

    if (name_size <= 19 && name.compare(0, name_size, "specific-activation", name_size) == 0)
      matches.push_back("specific-activation");

    if (name_size <= 5 && name.compare(0, name_size, "total", name_size) == 0)
      matches.push_back("total");

    if (name_size <= 7 && name.compare(0, name_size, "version", name_size) == 0)
      matches.push_back("version");

    if (name_size <= 4 && name.compare(0, name_size, "help", name_size) == 0)
      matches.push_back("help");

    if (matches.empty())
      {
        std::string error("unknown option '"); error += name; error += "'";
        throw clo::option_error(error);
      }

    if (matches.size() == 1)
      {
        return matches[0];
      }

    std::string error("the option name '"); error += name; error += "' is ambiguous";
    throw clo::option_error(error);
  }
} // end anonymous namespace