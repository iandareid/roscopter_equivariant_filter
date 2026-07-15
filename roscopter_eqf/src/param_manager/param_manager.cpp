#include "param_manager/param_manager.hpp"

#include <stdexcept>

namespace roscopter
{

ParamManager::ParamManager(rclcpp::Node * node)
: container_node_(node)
{
}

void ParamManager::declare_double(std::string param_name, const double value)
{
  params_[param_name] = value;
  container_node_->declare_parameter(param_name, value);
}

void ParamManager::declare_bool(std::string param_name, const bool value)
{
  params_[param_name] = value;
  container_node_->declare_parameter(param_name, value);
}

void ParamManager::declare_int(std::string param_name, const int64_t value)
{
  params_[param_name] = value;
  container_node_->declare_parameter(param_name, value);
}

void ParamManager::declare_string(std::string param_name, std::string value)
{
  params_[param_name] = value;
  container_node_->declare_parameter(param_name, value);
}

void ParamManager::set_double(std::string param_name, double value)
{
  if (params_.find(param_name) == params_.end()) {
    RCLCPP_ERROR_STREAM(container_node_->get_logger(), "Parameter not found: " << param_name);
    return;
  }
  params_[param_name] = value;
  container_node_->set_parameter(rclcpp::Parameter(param_name, value));
}

void ParamManager::set_bool(std::string param_name, bool value)
{
  if (params_.find(param_name) == params_.end()) {
    RCLCPP_ERROR_STREAM(container_node_->get_logger(), "Parameter not found: " << param_name);
    return;
  }
  params_[param_name] = value;
  container_node_->set_parameter(rclcpp::Parameter(param_name, value));
}

void ParamManager::set_int(std::string param_name, int64_t value)
{
  if (params_.find(param_name) == params_.end()) {
    RCLCPP_ERROR_STREAM(container_node_->get_logger(), "Parameter not found: " << param_name);
    return;
  }
  params_[param_name] = value;
  container_node_->set_parameter(rclcpp::Parameter(param_name, value));
}

void ParamManager::set_string(std::string param_name, std::string value)
{
  if (params_.find(param_name) == params_.end()) {
    RCLCPP_ERROR_STREAM(container_node_->get_logger(), "Parameter not found: " << param_name);
    return;
  }
  params_[param_name] = value;
  container_node_->set_parameter(rclcpp::Parameter(param_name, value));
}

double ParamManager::get_double(std::string param_name)
{
  return get_value<double>(param_name);
}

bool ParamManager::get_bool(std::string param_name)
{
  return get_value<bool>(param_name);
}

int64_t ParamManager::get_int(std::string param_name)
{
  return get_value<int64_t>(param_name);
}

std::string ParamManager::get_string(std::string param_name)
{
  return get_value<std::string>(param_name);
}

void ParamManager::set_parameters()
{
  for (auto & [key, value] : params_) {
    const rclcpp::Parameter parameter = container_node_->get_parameter(key);

    if (parameter.get_type() == rclcpp::ParameterType::PARAMETER_DOUBLE) {
      value = parameter.as_double();
    } else if (parameter.get_type() == rclcpp::ParameterType::PARAMETER_BOOL) {
      value = parameter.as_bool();
    } else if (parameter.get_type() == rclcpp::ParameterType::PARAMETER_INTEGER) {
      value = parameter.as_int();
    } else if (parameter.get_type() == rclcpp::ParameterType::PARAMETER_STRING) {
      value = parameter.as_string();
    } else {
      RCLCPP_ERROR_STREAM(
        container_node_->get_logger(),
        "Unable to set parameter: " << key);
    }
  }
}

bool ParamManager::set_parameters_callback(const std::vector<rclcpp::Parameter> & parameters)
{
  for (const auto & parameter : parameters) {
    const std::string name = parameter.get_name();
    if (params_.find(name) == params_.end()) {
      RCLCPP_ERROR_STREAM(
        container_node_->get_logger(),
        "Parameter is not declared in this node: " << name);
      return false;
    }

    if (parameter.get_type() == rclcpp::ParameterType::PARAMETER_DOUBLE) {
      params_[name] = parameter.as_double();
    } else if (parameter.get_type() == rclcpp::ParameterType::PARAMETER_BOOL) {
      params_[name] = parameter.as_bool();
    } else if (parameter.get_type() == rclcpp::ParameterType::PARAMETER_INTEGER) {
      params_[name] = parameter.as_int();
    } else if (parameter.get_type() == rclcpp::ParameterType::PARAMETER_STRING) {
      params_[name] = parameter.as_string();
    } else {
      RCLCPP_ERROR_STREAM(
        container_node_->get_logger(),
        "Unsupported parameter type for: " << name);
      return false;
    }
  }

  return true;
}

template<typename T>
T ParamManager::get_value(std::string param_name)
{
  const auto param = params_.find(param_name);
  if (param == params_.end()) {
    RCLCPP_ERROR_STREAM(container_node_->get_logger(), "Parameter not found: " << param_name);
    throw std::runtime_error("Parameter not found: " + param_name);
  }

  try {
    return std::get<T>(param->second);
  } catch (const std::bad_variant_access & exception) {
    RCLCPP_ERROR_STREAM(container_node_->get_logger(), "Parameter has wrong type: " << param_name);
    throw std::runtime_error(exception.what());
  }
}

template double ParamManager::get_value<double>(std::string);
template bool ParamManager::get_value<bool>(std::string);
template int64_t ParamManager::get_value<int64_t>(std::string);
template std::string ParamManager::get_value<std::string>(std::string);

}  // namespace roscopter
