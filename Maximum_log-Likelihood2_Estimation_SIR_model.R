{
  "cells": [
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "SDwZ7DSEkSEs"
      },
      "source": [
        "# Performing maximum likelihood estimation"
      ]
    },
    {
      "cell_type": "code",
      "source": [
        "install.packages(\"deSolve\")"
      ],
      "metadata": {
        "colab": {
          "base_uri": "https://localhost:8080/"
        },
        "id": "VnQhofALk0fM",
        "outputId": "72411e07-40bd-4980-fa12-10190bc6ad97"
      },
      "execution_count": null,
      "outputs": [
        {
          "output_type": "stream",
          "name": "stderr",
          "text": [
            "Installing package into ‘/usr/local/lib/R/site-library’\n",
            "(as ‘lib’ is unspecified)\n",
            "\n"
          ]
        }
      ]
    },
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {
        "id": "xSr_A1dBkSEx"
      },
      "outputs": [],
      "source": [
        "# set working directory >> this should be the directory that contains your code and data\n",
        "\n",
        "# load packages\n",
        "library(tidyverse)\n",
        "library(ggplot2)\n",
        "library(lubridate)\n",
        "library(scales)"
      ]
    },
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {
        "id": "Afj_qEG2kSE2"
      },
      "outputs": [],
      "source": [
        "# load data\n",
        "data <- read_csv(\"w.csv\", col_names = TRUE, show_col_types = FALSE)\n",
        "head(data)\n",
        "# keep only prevalence data\n",
        "data_cases <- select(data, c(\"Date_1\",  \"ACTIVE CASES\"))\n",
        "names(data_cases) <- c(\"date\", \"I\")   # rename columns\n",
        "\n",
        "# add a numeric data column\n",
        "#numeric date (counts of days from the starting date)\n",
        "d <- as.numeric(data_cases$date) \n",
        "d <- d - head(d,1) \n",
        "data_cases$time <- d\n",
        "\n",
        "# number of columns with missing values\n",
        "print(is.na(data_cases) %>% sum())\n",
        "\n",
        "# drop missing values\n",
        "data_cases <- data_cases %>% drop_na()\n",
        "\n",
        "# visualise data\n",
        "plot(data_cases$time, data_cases$I,\n",
        "    xlab = \"time (days)\", ylab = \"number of active cases\", type='line')"
      ]
    },
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "sxqqv8kDkSE3"
      },
      "source": [
        "You can see that the data is so incomplete. Moreover, a SIR model cannot capture the occurence of multiple infection waves. So, probably, it is better to fit to data for the first wave only."
      ]
    },
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {
        "id": "gCgh0kvTkSE3"
      },
      "outputs": [],
      "source": [
        "# data for first wave, assuming it ends at t = 125 days from the starting date\n",
        "dat_wave_1 <- data_cases[data_cases$time < 125 & data_cases$I > 0,]\n",
        "dat_wave_1$I <- as.numeric(dat_wave_1$I)   # convert I to numeric format\n",
        "\n",
        "head(dat_wave_1)"
      ]
    },
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {
        "id": "eEqXtxFGkSE5"
      },
      "outputs": [],
      "source": [
        "# log-likelihood function\n",
        "SIR_LLH <- function(params,           # model parameters\n",
        "                    model_func,       # model function\n",
        "                    init_cond,        # initial condition\n",
        "                    time_step = .1,   # integration step\n",
        "                    r = 1,            # case reporting rate\n",
        "                    data) {           # data\n",
        "    \n",
        "    # Calculate model output\n",
        "    # load packages\n",
        "    require(deSolve)\n",
        "    \n",
        "    # integration time points\n",
        "    times <- seq(0, max(data$time), by = time_step)  \n",
        "    \n",
        "    # simulate the model with given parameter values\n",
        "    output <- as.data.frame(ode(y = init_cond,\n",
        "                               times = times,\n",
        "                               func = model_func,\n",
        "                               parms = params))\n",
        "    \n",
        "    # Calculate log-likelihood (LLH) of model fit\n",
        "    data <- na.omit(data)\n",
        "    \n",
        "    LLH <- sum(dpois(data$I, r * output$I[output$time %in% data$time], log = TRUE))\n",
        "\n",
        "    return(LLH)\n",
        "\n",
        "}\n",
        "\n",
        "#----------------------#\n",
        "# sir model function   #\n",
        "#----------------------#\n",
        "\n",
        "sir_model <- function(time, state, params){\n",
        "    with(as.list(c(state, params)), {\n",
        "        N <- S+I+R           # population size\n",
        "        lambda <- beta * I / N  # force of infection\n",
        "        \n",
        "        # model equations\n",
        "        dS <- - lambda * S \n",
        "        dI <- lambda * S - gamma * I\n",
        "        dR <- gamma * I\n",
        "            \n",
        "        return(list(c(dS, dI, dR)))\n",
        "    })\n",
        "}\n"
      ]
    },
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {
        "id": "eaq2Zu_jkSE_"
      },
      "outputs": [],
      "source": [
        "# OPTIMISATION:\n",
        "## assume the total population is N = 60 million and I(0) = 1\n",
        "N <- 60041996\n",
        "initial_state_values <- c(S = N-1, I = 1, R = 0)\n",
        "\n",
        "optimised <- optim(par = c(beta = 0.125, gamma = 0.071), # initialvalues for gamma and beta\n",
        "                    fn = SIR_LLH,\n",
        "                    dat = dat_wave_1,\n",
        "                    model_func = sir_model,\n",
        "                    init_cond = initial_state_values,\n",
        "                    r = 1,\n",
        "                    control = list(fnscale=-1))  # tells optim() to look for the maximum number instead of the minimum (the default)\n",
        "                    \n",
        "optimised$par"
      ]
    },
    {
      "cell_type": "markdown",
      "metadata": {
        "id": "66Pql-tXkSFA"
      },
      "source": [
        "Finally, since we have the full dataset of the number of infected people, confirm that these parameter values indeed produce a good visual fit to the real data by plotting the model simulation alongside the data, in the cell below."
      ]
    },
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {
        "id": "_njGasdSkSFA"
      },
      "outputs": [],
      "source": [
        "# compute best fit solution\n",
        "times <- seq(0, max(dat_wave_1$time), by = .5)\n",
        "\n",
        "best_fit <- as.data.frame(ode(y = initial_state_values,\n",
        "                              times = times,\n",
        "                              func = sir_model,\n",
        "                              parms = optimised$par))\n",
        "            \n",
        "require(ggplot2)\n",
        "options(repr.plot.width = 12, repr.plot.height = 8)\n",
        "ggplot() + \n",
        "geom_point(data = dat_wave_1, aes(x = time, y = I, color = 'data')) +\n",
        "geom_line(data = best_fit, aes(x = time, y = I,  color = \"SIR model\")) +\n",
        "xlab('time (days)') +\n",
        "ylab('number of cases') +\n",
        "labs(title = paste(\"Likelihood fit of an SIR model \"))"
      ]
    },
    {
      "cell_type": "code",
      "execution_count": null,
      "metadata": {
        "id": "NnsOURO4kSFC"
      },
      "outputs": [],
      "source": [
        ""
      ]
    }
  ],
  "metadata": {
    "kernelspec": {
      "display_name": "R",
      "language": "R",
      "name": "ir"
    },
    "language_info": {
      "codemirror_mode": "r",
      "file_extension": ".r",
      "mimetype": "text/x-r-source",
      "name": "R",
      "pygments_lexer": "r",
      "version": "4.2.0"
    },
    "colab": {
      "name": "Maximum_Likelihood2_Estimation_SIR_model.ipynb",
      "provenance": []
    }
  },
  "nbformat": 4,
  "nbformat_minor": 0
}